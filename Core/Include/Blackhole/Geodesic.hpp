/**
 * @file geodesic.hpp
 * @brief Geodesic Integrator for Null and Timelike Geodesics
 * 
 * Implements high-order Runge-Kutta integrators (RK4, RKF45) with adaptive
 * step sizing for ray tracing in curved spacetime.
 * 
 * @author BlackHole-Sim Team
 * @license MIT
 */

#ifndef BLACKHOLE_GEODESIC_HPP
#define BLACKHOLE_GEODESIC_HPP

#include "metric.hpp"
#include <vector>
#include <functional>
#include <optional>
#include <algorithm>

namespace blackhole {

/**
 * @brief Result of geodesic integration
 */
template<typename T = double>
struct GeodesicResult {
    std::vector<Vec4<T>> positions;     // Position along geodesic
    std::vector<Vec4<T>> momenta;       // 4-momentum along geodesic
    std::vector<T> affine_params;       // Affine parameter values
    
    bool hit_horizon = false;           // Did ray cross event horizon?
    bool escaped = false;               // Did ray escape to infinity?
    int num_steps = 0;                  // Number of integration steps
    
    // Final state
    Vec4<T> final_position;
    Vec4<T> final_momentum;
    T final_affine_param = 0;
};

/**
 * @brief State vector for geodesic integration (position + momentum)
 */
template<typename T = double>
struct GeodesicState {
    Vec4<T> x;  // Position (t, r, θ, φ)
    Vec4<T> p;  // 4-momentum (p_t, p_r, p_θ, p_φ)
    
    GeodesicState operator+(const GeodesicState& other) const {
        return {
            Vec4<T>(x.t + other.x.t, x.r + other.x.r, x.theta + other.x.theta, x.phi + other.x.phi),
            Vec4<T>(p.t + other.p.t, p.r + other.p.r, p.theta + other.p.theta, p.phi + other.p.phi)
        };
    }
    
    GeodesicState operator*(T scalar) const {
        return {
            Vec4<T>(x.t * scalar, x.r * scalar, x.theta * scalar, x.phi * scalar),
            Vec4<T>(p.t * scalar, p.r * scalar, p.theta * scalar, p.phi * scalar)
        };
    }
};

/**
 * @brief Configuration for geodesic integration
 */
template<typename T = double>
struct IntegratorConfig {
    T max_affine_param = 1000.0;        // Maximum affine parameter
    T initial_step = 0.1;               // Initial step size
    T min_step = 1e-8;                  // Minimum step size
    T max_step = 10.0;                  // Maximum step size
    T tolerance = 1e-8;                 // Error tolerance for adaptive stepping
    int max_steps = 100000;             // Maximum number of steps
    T horizon_threshold = 1.01;         // Stop when r < threshold * r_horizon
    T escape_radius = 1000.0;           // Consider escaped beyond this radius
    bool store_trajectory = false;      // Whether to store full trajectory
    bool backward = true;               // Integrate backward (for ray tracing)
};

/**
 * @brief Geodesic integrator using the geodesic equation
 * 
 * Integrates: dx^μ/dλ = p^μ
 *             dp^μ/dλ = -Γ^μ_νρ p^ν p^ρ
 */
template<typename T = double>
class GeodesicIntegrator {
private:
    const Metric<T>& metric_;
    IntegratorConfig<T> config_;
    
    /**
     * @brief Compute right-hand side of geodesic equations
     */
    GeodesicState<T> geodesic_rhs(const GeodesicState<T>& state) const {
        GeodesicState<T> dstate;
        
        // dx^μ/dλ = g^μν p_ν (raise index of momentum)
        auto g = metric_.compute(state.x.r, state.x.theta);
        
        for (int mu = 0; mu < 4; ++mu) {
            dstate.x[mu] = 0;
            for (int nu = 0; nu < 4; ++nu) {
                dstate.x[mu] += g.g_inv[mu][nu] * state.p[nu];
            }
        }
        
        // dp_μ/dλ = -Γ^ν_μρ p_ν p^ρ = -(1/2) ∂g_νρ/∂x^μ p^ν p^ρ
        // Using Christoffel symbols: dp^μ/dλ = -Γ^μ_νρ p^ν p^ρ
        auto Gamma = metric_.christoffel(state.x.r, state.x.theta);
        
        // Get contravariant momentum
        Vec4<T> p_up;
        for (int mu = 0; mu < 4; ++mu) {
            p_up[mu] = 0;
            for (int nu = 0; nu < 4; ++nu) {
                p_up[mu] += g.g_inv[mu][nu] * state.p[nu];
            }
        }
        
        for (int mu = 0; mu < 4; ++mu) {
            dstate.p[mu] = 0;
            for (int nu = 0; nu < 4; ++nu) {
                for (int rho = 0; rho < 4; ++rho) {
                    // dp_μ/dλ = Γ^ν_μρ p_ν p^ρ (note sign for covariant momentum)
                    dstate.p[mu] += Gamma[nu][mu][rho] * state.p[nu] * p_up[rho];
                }
            }
        }
        
        return dstate;
    }
    
    /**
     * @brief Single RK4 step
     */
    GeodesicState<T> rk4_step(const GeodesicState<T>& state, T h) const {
        auto k1 = geodesic_rhs(state);
        auto k2 = geodesic_rhs(state + k1 * (h / 2));
        auto k3 = geodesic_rhs(state + k2 * (h / 2));
        auto k4 = geodesic_rhs(state + k3 * h);
        
        return state + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0);
    }
    
    /**
     * @brief RK45 (Dormand-Prince) step with error estimation
     * @return {new_state, error_estimate}
     */
    std::pair<GeodesicState<T>, T> rkf45_step(const GeodesicState<T>& state, T h) const {
        // Dormand-Prince coefficients
        constexpr T a2 = 1.0/5.0, a3 = 3.0/10.0, a4 = 4.0/5.0, a5 = 8.0/9.0;
        constexpr T b21 = 1.0/5.0;
        constexpr T b31 = 3.0/40.0, b32 = 9.0/40.0;
        constexpr T b41 = 44.0/45.0, b42 = -56.0/15.0, b43 = 32.0/9.0;
        constexpr T b51 = 19372.0/6561.0, b52 = -25360.0/2187.0, b53 = 64448.0/6561.0, b54 = -212.0/729.0;
        constexpr T b61 = 9017.0/3168.0, b62 = -355.0/33.0, b63 = 46732.0/5247.0, b64 = 49.0/176.0, b65 = -5103.0/18656.0;
        
        constexpr T c1 = 35.0/384.0, c3 = 500.0/1113.0, c4 = 125.0/192.0, c5 = -2187.0/6784.0, c6 = 11.0/84.0;
        constexpr T d1 = 5179.0/57600.0, d3 = 7571.0/16695.0, d4 = 393.0/640.0, d5 = -92097.0/339200.0, d6 = 187.0/2100.0, d7 = 1.0/40.0;
        
        auto k1 = geodesic_rhs(state);
        auto k2 = geodesic_rhs(state + k1 * (h * b21));
        auto k3 = geodesic_rhs(state + k1 * (h * b31) + k2 * (h * b32));
        auto k4 = geodesic_rhs(state + k1 * (h * b41) + k2 * (h * b42) + k3 * (h * b43));
        auto k5 = geodesic_rhs(state + k1 * (h * b51) + k2 * (h * b52) + k3 * (h * b53) + k4 * (h * b54));
        auto k6 = geodesic_rhs(state + k1 * (h * b61) + k2 * (h * b62) + k3 * (h * b63) + k4 * (h * b64) + k5 * (h * b65));
        
        // 5th order solution
        auto new_state = state + (k1 * c1 + k3 * c3 + k4 * c4 + k5 * c5 + k6 * c6) * h;
        
        // 4th order solution for error estimate
        auto k7 = geodesic_rhs(new_state);
        auto err_state = (k1 * (c1 - d1) + k3 * (c3 - d3) + k4 * (c4 - d4) + 
                         k5 * (c5 - d5) + k6 * (c6 - d6) + k7 * (-d7)) * h;
        
        // Error estimate (use r and theta components as most sensitive)
        T error = std::sqrt(err_state.x.r * err_state.x.r + 
                           err_state.x.theta * err_state.x.theta +
                           err_state.p.r * err_state.p.r +
                           err_state.p.theta * err_state.p.theta);
        
        return {new_state, error};
    }
    
public:
    GeodesicIntegrator(const Metric<T>& metric, const IntegratorConfig<T>& config = {})
        : metric_(metric), config_(config) {}
    
    /**
     * @brief Initialize null geodesic from camera pixel
     * 
     * @param r_obs Observer radial position
     * @param theta_obs Observer polar angle
     * @param alpha Impact parameter (horizontal angle)
     * @param beta Impact parameter (vertical angle)
     * @return Initial state for geodesic
     */
    GeodesicState<T> init_null_geodesic(T r_obs, T theta_obs, T alpha, T beta) const {
        GeodesicState<T> state;
        
        // Observer position
        state.x.t = 0;
        state.x.r = r_obs;
        state.x.theta = theta_obs;
        state.x.phi = 0;
        
        // Compute initial 4-momentum for null geodesic
        auto g = metric_.compute(r_obs, theta_obs);
        
        T sin_theta = std::sin(theta_obs);
        T cos_theta = std::cos(theta_obs);
        
        // Direction in local frame
        T k_r = -std::cos(beta) * std::cos(alpha);
        T k_theta = std::sin(beta);
        T k_phi = std::cos(beta) * std::sin(alpha);
        
        // Transform to coordinate basis
        state.p.r = k_r * std::sqrt(g.g[1][1]);
        state.p.theta = k_theta * std::sqrt(g.g[2][2]);
        state.p.phi = k_phi * std::sqrt(g.g[3][3]) / sin_theta;
        
        // For null geodesic: g^μν p_μ p_ν = 0
        // Solve for p_t
        T A = g.g_inv[0][0];
        T B = 2.0 * g.g_inv[0][3] * state.p.phi;
        T C = g.g_inv[1][1] * state.p.r * state.p.r + 
              g.g_inv[2][2] * state.p.theta * state.p.theta +
              g.g_inv[3][3] * state.p.phi * state.p.phi;
        
        // p_t = (-B ± sqrt(B² - 4AC)) / 2A
        T discriminant = B * B - 4.0 * A * C;
        if (discriminant < 0) {
            discriminant = 0; // Numerical safety
        }
        
        // Choose sign for ingoing ray (negative energy at infinity)
        state.p.t = (-B - std::sqrt(discriminant)) / (2.0 * A);
        
        return state;
    }
    
    /**
     * @brief Integrate geodesic from initial state
     */
    GeodesicResult<T> integrate(GeodesicState<T> state) const {
        GeodesicResult<T> result;
        
        T lambda = 0;
        T h = config_.initial_step;
        if (config_.backward) h = -h;
        
        T r_horizon = metric_.horizon_radius() * config_.horizon_threshold;
        
        if (config_.store_trajectory) {
            result.positions.push_back(state.x);
            result.momenta.push_back(state.p);
            result.affine_params.push_back(lambda);
        }
        
        while (result.num_steps < config_.max_steps) {
            // Check termination conditions
            if (state.x.r < r_horizon) {
                result.hit_horizon = true;
                break;
            }
            
            if (state.x.r > config_.escape_radius) {
                result.escaped = true;
                break;
            }
            
            if (std::abs(lambda) > config_.max_affine_param) {
                break;
            }
            
            // Keep theta in bounds
            if (state.x.theta < 0.01) state.x.theta = 0.01;
            if (state.x.theta > constants::pi - 0.01) state.x.theta = constants::pi - 0.01;
            
            // Adaptive RKF45 step
            auto [new_state, error] = rkf45_step(state, h);
            
            // Adjust step size
            if (error > config_.tolerance && std::abs(h) > config_.min_step) {
                h *= 0.5;
                continue;
            }
            
            // Accept step
            state = new_state;
            lambda += h;
            result.num_steps++;
            
            if (config_.store_trajectory) {
                result.positions.push_back(state.x);
                result.momenta.push_back(state.p);
                result.affine_params.push_back(lambda);
            }
            
            // Grow step size if error is small
            if (error < config_.tolerance * 0.1 && std::abs(h) < config_.max_step) {
                h *= 2.0;
            }
            
            // Limit step size near horizon
            T r_factor = (state.x.r - r_horizon) / metric_.mass();
            if (r_factor < 10.0 && r_factor > 0) {
                h = std::copysign(std::min(std::abs(h), r_factor * 0.1), h);
            }
        }
        
        result.final_position = state.x;
        result.final_momentum = state.p;
        result.final_affine_param = lambda;
        
        return result;
    }
    
    /**
     * @brief Trace ray from observer to source (backward ray tracing)
     */
    GeodesicResult<T> trace_ray(T r_obs, T theta_obs, T alpha, T beta) const {
        auto state = init_null_geodesic(r_obs, theta_obs, alpha, beta);
        return integrate(state);
    }
    
    /**
     * @brief Check if a ray hits the accretion disk
     * @param result Geodesic integration result
     * @param r_in Inner edge of disk
     * @param r_out Outer edge of disk
     * @return Optional<{r, phi}> if ray crosses disk equatorial plane
     */
    std::optional<std::pair<T, T>> disk_intersection(
        const GeodesicResult<T>& result, T r_in, T r_out) const {
        
        if (result.positions.size() < 2) return std::nullopt;
        
        // Look for equatorial plane crossing
        for (size_t i = 1; i < result.positions.size(); ++i) {
            T theta_prev = result.positions[i-1].theta;
            T theta_curr = result.positions[i].theta;
            
            T theta_eq = constants::pi / 2.0;
            
            // Check if we crossed the equatorial plane
            if ((theta_prev - theta_eq) * (theta_curr - theta_eq) < 0) {
                // Linear interpolation to find crossing point
                T t = (theta_eq - theta_prev) / (theta_curr - theta_prev);
                T r = result.positions[i-1].r + t * (result.positions[i].r - result.positions[i-1].r);
                T phi = result.positions[i-1].phi + t * (result.positions[i].phi - result.positions[i-1].phi);
                
                if (r >= r_in && r <= r_out) {
                    return std::make_pair(r, phi);
                }
            }
        }
        
        return std::nullopt;
    }
};

/**
 * @brief Camera model for generating image rays
 */
template<typename T = double>
class Camera {
private:
    T r_obs_;           // Observer radial coordinate
    T theta_obs_;       // Observer polar angle
    T fov_;             // Field of view (radians)
    int width_;         // Image width (pixels)
    int height_;        // Image height (pixels)
    
public:
    Camera(T r_obs, T theta_obs, T fov_deg, int width, int height)
        : r_obs_(r_obs), theta_obs_(theta_obs), 
          fov_(fov_deg * constants::pi / 180.0),
          width_(width), height_(height) {}
    
    /**
     * @brief Get impact parameters for a pixel
     * @param i Pixel column (0 to width-1)
     * @param j Pixel row (0 to height-1)
     * @return {alpha, beta} impact parameters
     */
    std::pair<T, T> pixel_to_angles(int i, int j) const {
        T aspect = static_cast<T>(width_) / height_;
        
        T x = (2.0 * i / width_ - 1.0) * std::tan(fov_ / 2.0) * aspect;
        T y = (1.0 - 2.0 * j / height_) * std::tan(fov_ / 2.0);
        
        T alpha = std::atan(x);
        T beta = std::atan(y / std::sqrt(1 + x * x));
        
        return {alpha, beta};
    }
    
    T r_observer() const { return r_obs_; }
    T theta_observer() const { return theta_obs_; }
    int width() const { return width_; }
    int height() const { return height_; }
};

/**
 * @brief Ray tracer for generating black hole images
 */
template<typename T = double>
class RayTracer {
private:
    const Metric<T>& metric_;
    const Camera<T>& camera_;
    IntegratorConfig<T> config_;
    
public:
    RayTracer(const Metric<T>& metric, const Camera<T>& camera,
              const IntegratorConfig<T>& config = {})
        : metric_(metric), camera_(camera), config_(config) {}
    
    /**
     * @brief Trace single pixel
     */
    GeodesicResult<T> trace_pixel(int i, int j) const {
        GeodesicIntegrator<T> integrator(metric_, config_);
        auto [alpha, beta] = camera_.pixel_to_angles(i, j);
        return integrator.trace_ray(camera_.r_observer(), camera_.theta_observer(), alpha, beta);
    }
    
    /**
     * @brief Generate full image (intensity values)
     * @param disk_inner Inner radius of accretion disk
     * @param disk_outer Outer radius of accretion disk
     * @return 2D vector of intensity values
     */
    std::vector<std::vector<T>> render(T disk_inner, T disk_outer) const {
        int w = camera_.width();
        int h = camera_.height();
        
        std::vector<std::vector<T>> image(h, std::vector<T>(w, 0.0));
        GeodesicIntegrator<T> integrator(metric_, config_);
        
        #pragma omp parallel for collapse(2) schedule(dynamic)
        for (int j = 0; j < h; ++j) {
            for (int i = 0; i < w; ++i) {
                auto [alpha, beta] = camera_.pixel_to_angles(i, j);
                auto result = integrator.trace_ray(
                    camera_.r_observer(), camera_.theta_observer(), alpha, beta);
                
                // Check for disk intersection
                auto hit = integrator.disk_intersection(result, disk_inner, disk_outer);
                
                if (hit) {
                    T r = hit->first;
                    T phi = hit->second;
                    
                    // Simple disk emission model (temperature ~ r^(-3/4))
                    T intensity = std::pow(disk_inner / r, 3.0);
                    
                    // Add doppler beaming approximation
                    T v_orbit = std::sqrt(1.0 / r);
                    T doppler = 1.0 + v_orbit * std::sin(phi) * std::sin(camera_.theta_observer());
                    intensity *= std::pow(std::abs(doppler), 3.0);
                    
                    image[j][i] = intensity;
                } else if (result.hit_horizon) {
                    image[j][i] = 0.0;  // Pure black for horizon
                } else if (result.escaped) {
                    // Background (could add star field)
                    image[j][i] = 0.01;
                }
            }
        }
        
        return image;
    }
};

} // namespace blackhole

#endif // BLACKHOLE_GEODESIC_HPP