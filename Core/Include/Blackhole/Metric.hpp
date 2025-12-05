/**
 * @file metric.hpp
 * @brief General Relativistic Metric Tensor Classes
 * 
 * Implements Schwarzschild and Kerr metrics with full tensor algebra,
 * Christoffel symbols, and coordinate transformations.
 * 
 * @author BlackHole-Sim Team
 * @license MIT
 */

#ifndef BLACKHOLE_METRIC_HPP
#define BLACKHOLE_METRIC_HPP

#include <cmath>
#include <array>
#include <stdexcept>
#include <memory>

namespace blackhole {

// Physical constants in CGS units
namespace constants {
    constexpr double c = 2.99792458e10;      // Speed of light (cm/s)
    constexpr double G = 6.67430e-8;          // Gravitational constant (cm³/g/s²)
    constexpr double M_sun = 1.98892e33;      // Solar mass (g)
    constexpr double pi = 3.14159265358979323846;
}

/**
 * @brief 4-vector in spacetime
 */
template<typename T = double>
struct Vec4 {
    T t, r, theta, phi;
    
    Vec4() : t(0), r(0), theta(0), phi(0) {}
    Vec4(T t_, T r_, T theta_, T phi_) : t(t_), r(r_), theta(theta_), phi(phi_) {}
    
    T& operator[](int i) {
        switch(i) {
            case 0: return t;
            case 1: return r;
            case 2: return theta;
            case 3: return phi;
            default: throw std::out_of_range("Vec4 index out of range");
        }
    }
    
    const T& operator[](int i) const {
        switch(i) {
            case 0: return t;
            case 1: return r;
            case 2: return theta;
            case 3: return phi;
            default: throw std::out_of_range("Vec4 index out of range");
        }
    }
    
    Vec4 operator+(const Vec4& v) const {
        return Vec4(t + v.t, r + v.r, theta + v.theta, phi + v.phi);
    }
    
    Vec4 operator*(T scalar) const {
        return Vec4(t * scalar, r * scalar, theta * scalar, phi * scalar);
    }
};

/**
 * @brief 4x4 symmetric metric tensor
 */
template<typename T = double>
class MetricTensor {
public:
    std::array<std::array<T, 4>, 4> g;      // Covariant metric g_μν
    std::array<std::array<T, 4>, 4> g_inv;  // Contravariant metric g^μν
    
    MetricTensor() {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                g[i][j] = (i == j) ? 1.0 : 0.0;
                g_inv[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
    }
    
    T& operator()(int mu, int nu) { return g[mu][nu]; }
    const T& operator()(int mu, int nu) const { return g[mu][nu]; }
    
    /**
     * @brief Compute determinant of the metric tensor
     */
    T determinant() const {
        // For diagonal + g_tφ metrics (Kerr in BL coords)
        T det = g[0][0] * g[1][1] * g[2][2] * g[3][3];
        det -= g[0][3] * g[0][3] * g[1][1] * g[2][2];
        return det;
    }
    
    /**
     * @brief Lower an index (contravariant → covariant)
     */
    Vec4<T> lower(const Vec4<T>& v) const {
        Vec4<T> result;
        for (int mu = 0; mu < 4; ++mu) {
            result[mu] = 0;
            for (int nu = 0; nu < 4; ++nu) {
                result[mu] += g[mu][nu] * v[nu];
            }
        }
        return result;
    }
    
    /**
     * @brief Raise an index (covariant → contravariant)
     */
    Vec4<T> raise(const Vec4<T>& v) const {
        Vec4<T> result;
        for (int mu = 0; mu < 4; ++mu) {
            result[mu] = 0;
            for (int nu = 0; nu < 4; ++nu) {
                result[mu] += g_inv[mu][nu] * v[nu];
            }
        }
        return result;
    }
    
    /**
     * @brief Compute inner product of two 4-vectors
     */
    T inner(const Vec4<T>& u, const Vec4<T>& v) const {
        T result = 0;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                result += g[mu][nu] * u[mu] * v[nu];
            }
        }
        return result;
    }
};

/**
 * @brief Abstract base class for spacetime metrics
 */
template<typename T = double>
class Metric {
public:
    virtual ~Metric() = default;
    
    /**
     * @brief Compute metric tensor at given coordinates
     * @param r Radial coordinate
     * @param theta Polar angle
     * @return MetricTensor at (r, θ)
     */
    virtual MetricTensor<T> compute(T r, T theta) const = 0;
    
    /**
     * @brief Compute Christoffel symbols Γ^α_μν at given coordinates
     * @param r Radial coordinate
     * @param theta Polar angle
     * @return 4x4x4 array of Christoffel symbols
     */
    virtual std::array<std::array<std::array<T, 4>, 4>, 4> 
    christoffel(T r, T theta) const = 0;
    
    /**
     * @brief Get mass parameter
     */
    virtual T mass() const = 0;
    
    /**
     * @brief Get spin parameter
     */
    virtual T spin() const = 0;
    
    /**
     * @brief Compute event horizon radius
     */
    virtual T horizon_radius() const = 0;
    
    /**
     * @brief Compute ISCO radius (prograde)
     */
    virtual T isco_radius() const = 0;
    
    /**
     * @brief Compute photon sphere radius (equatorial, prograde)
     */
    virtual T photon_sphere_radius() const = 0;
};

/**
 * @brief Schwarzschild metric (non-rotating black hole)
 * 
 * ds² = -(1 - 2M/r)dt² + (1 - 2M/r)⁻¹dr² + r²dθ² + r²sin²θ dφ²
 */
template<typename T = double>
class SchwarzschildMetric : public Metric<T> {
private:
    T M_;  // Mass parameter (in geometric units where G = c = 1)
    
public:
    explicit SchwarzschildMetric(T mass = 1.0) : M_(mass) {
        if (mass <= 0) {
            throw std::invalid_argument("Mass must be positive");
        }
    }
    
    T mass() const override { return M_; }
    T spin() const override { return 0.0; }
    
    T horizon_radius() const override { 
        return 2.0 * M_; 
    }
    
    T isco_radius() const override { 
        return 6.0 * M_; 
    }
    
    T photon_sphere_radius() const override { 
        return 3.0 * M_; 
    }
    
    MetricTensor<T> compute(T r, T theta) const override {
        MetricTensor<T> metric;
        
        if (r <= 0) {
            throw std::domain_error("Radial coordinate must be positive");
        }
        
        T f = 1.0 - 2.0 * M_ / r;
        T sin_theta = std::sin(theta);
        T sin2_theta = sin_theta * sin_theta;
        
        // Covariant metric g_μν
        metric.g[0][0] = -f;           // g_tt
        metric.g[1][1] = 1.0 / f;      // g_rr
        metric.g[2][2] = r * r;        // g_θθ
        metric.g[3][3] = r * r * sin2_theta;  // g_φφ
        
        // Contravariant metric g^μν
        metric.g_inv[0][0] = -1.0 / f;
        metric.g_inv[1][1] = f;
        metric.g_inv[2][2] = 1.0 / (r * r);
        metric.g_inv[3][3] = 1.0 / (r * r * sin2_theta);
        
        return metric;
    }
    
    std::array<std::array<std::array<T, 4>, 4>, 4> 
    christoffel(T r, T theta) const override {
        std::array<std::array<std::array<T, 4>, 4>, 4> Gamma = {};
        
        T f = 1.0 - 2.0 * M_ / r;
        T r2 = r * r;
        T sin_theta = std::sin(theta);
        T cos_theta = std::cos(theta);
        T sin2_theta = sin_theta * sin_theta;
        
        // Non-zero Christoffel symbols
        // Γ^t_tr = Γ^t_rt
        Gamma[0][0][1] = M_ / (r2 * f);
        Gamma[0][1][0] = Gamma[0][0][1];
        
        // Γ^r_tt
        Gamma[1][0][0] = M_ * f / r2;
        
        // Γ^r_rr
        Gamma[1][1][1] = -M_ / (r2 * f);
        
        // Γ^r_θθ
        Gamma[1][2][2] = -(r - 2.0 * M_);
        
        // Γ^r_φφ
        Gamma[1][3][3] = -(r - 2.0 * M_) * sin2_theta;
        
        // Γ^θ_rθ = Γ^θ_θr
        Gamma[2][1][2] = 1.0 / r;
        Gamma[2][2][1] = Gamma[2][1][2];
        
        // Γ^θ_φφ
        Gamma[2][3][3] = -sin_theta * cos_theta;
        
        // Γ^φ_rφ = Γ^φ_φr
        Gamma[3][1][3] = 1.0 / r;
        Gamma[3][3][1] = Gamma[3][1][3];
        
        // Γ^φ_θφ = Γ^φ_φθ
        Gamma[3][2][3] = cos_theta / sin_theta;
        Gamma[3][3][2] = Gamma[3][2][3];
        
        return Gamma;
    }
};

/**
 * @brief Kerr metric (rotating black hole) in Boyer-Lindquist coordinates
 * 
 * ds² = -(1 - 2Mr/Σ)dt² - (4Mar sin²θ/Σ)dtdφ 
 *       + (Σ/Δ)dr² + Σdθ² + (A sin²θ/Σ)dφ²
 * 
 * where:
 *   Σ = r² + a²cos²θ
 *   Δ = r² - 2Mr + a²
 *   A = (r² + a²)² - a²Δsin²θ
 */
template<typename T = double>
class KerrMetric : public Metric<T> {
private:
    T M_;  // Mass parameter
    T a_;  // Spin parameter (|a| < M for black hole)
    
    // Helper functions
    T Sigma(T r, T theta) const {
        T cos_theta = std::cos(theta);
        return r * r + a_ * a_ * cos_theta * cos_theta;
    }
    
    T Delta(T r) const {
        return r * r - 2.0 * M_ * r + a_ * a_;
    }
    
    T A_func(T r, T theta) const {
        T r2_a2 = r * r + a_ * a_;
        T sin_theta = std::sin(theta);
        return r2_a2 * r2_a2 - a_ * a_ * Delta(r) * sin_theta * sin_theta;
    }
    
public:
    KerrMetric(T mass = 1.0, T spin = 0.0) : M_(mass), a_(spin) {
        if (mass <= 0) {
            throw std::invalid_argument("Mass must be positive");
        }
        if (std::abs(spin) >= mass) {
            throw std::invalid_argument("Spin magnitude must be less than mass for black hole");
        }
    }
    
    T mass() const override { return M_; }
    T spin() const override { return a_; }
    
    /**
     * @brief Outer event horizon radius r_+
     */
    T horizon_radius() const override {
        return M_ + std::sqrt(M_ * M_ - a_ * a_);
    }
    
    /**
     * @brief Inner (Cauchy) horizon radius r_-
     */
    T inner_horizon_radius() const {
        return M_ - std::sqrt(M_ * M_ - a_ * a_);
    }
    
    /**
     * @brief ISCO radius for prograde orbits
     */
    T isco_radius() const override {
        T a_star = a_ / M_;
        T Z1 = 1.0 + std::cbrt(1.0 - a_star * a_star) * 
               (std::cbrt(1.0 + a_star) + std::cbrt(1.0 - a_star));
        T Z2 = std::sqrt(3.0 * a_star * a_star + Z1 * Z1);
        return M_ * (3.0 + Z2 - std::sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));
    }
    
    /**
     * @brief ISCO radius for retrograde orbits
     */
    T isco_radius_retrograde() const {
        T a_star = a_ / M_;
        T Z1 = 1.0 + std::cbrt(1.0 - a_star * a_star) * 
               (std::cbrt(1.0 + a_star) + std::cbrt(1.0 - a_star));
        T Z2 = std::sqrt(3.0 * a_star * a_star + Z1 * Z1);
        return M_ * (3.0 + Z2 + std::sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));
    }
    
    /**
     * @brief Photon sphere radius (equatorial, prograde)
     */
    T photon_sphere_radius() const override {
        T a_star = a_ / M_;
        return 2.0 * M_ * (1.0 + std::cos(2.0/3.0 * std::acos(-a_star)));
    }
    
    /**
     * @brief Ergosphere outer boundary at given theta
     */
    T ergosphere_radius(T theta) const {
        T cos_theta = std::cos(theta);
        return M_ + std::sqrt(M_ * M_ - a_ * a_ * cos_theta * cos_theta);
    }
    
    /**
     * @brief Angular velocity at the horizon
     */
    T horizon_angular_velocity() const {
        T r_plus = horizon_radius();
        return a_ / (2.0 * M_ * r_plus);
    }
    
    MetricTensor<T> compute(T r, T theta) const override {
        MetricTensor<T> metric;
        
        if (r <= 0) {
            throw std::domain_error("Radial coordinate must be positive");
        }
        
        T sigma = Sigma(r, theta);
        T delta = Delta(r);
        T A = A_func(r, theta);
        T sin_theta = std::sin(theta);
        T sin2_theta = sin_theta * sin_theta;
        T cos_theta = std::cos(theta);
        
        // Covariant metric g_μν
        metric.g[0][0] = -(1.0 - 2.0 * M_ * r / sigma);
        metric.g[0][3] = -2.0 * M_ * a_ * r * sin2_theta / sigma;
        metric.g[3][0] = metric.g[0][3];
        metric.g[1][1] = sigma / delta;
        metric.g[2][2] = sigma;
        metric.g[3][3] = A * sin2_theta / sigma;
        
        // Contravariant metric g^μν
        T det_g_tphi = metric.g[0][0] * metric.g[3][3] - metric.g[0][3] * metric.g[0][3];
        
        metric.g_inv[0][0] = -metric.g[3][3] / det_g_tphi;
        metric.g_inv[0][3] = metric.g[0][3] / det_g_tphi;
        metric.g_inv[3][0] = metric.g_inv[0][3];
        metric.g_inv[1][1] = delta / sigma;
        metric.g_inv[2][2] = 1.0 / sigma;
        metric.g_inv[3][3] = -metric.g[0][0] / det_g_tphi;
        
        return metric;
    }
    
    std::array<std::array<std::array<T, 4>, 4>, 4> 
    christoffel(T r, T theta) const override {
        std::array<std::array<std::array<T, 4>, 4>, 4> Gamma = {};
        
        T sigma = Sigma(r, theta);
        T delta = Delta(r);
        T sin_theta = std::sin(theta);
        T cos_theta = std::cos(theta);
        T sin2_theta = sin_theta * sin_theta;
        T cos2_theta = cos_theta * cos_theta;
        
        T r2 = r * r;
        T a2 = a_ * a_;
        T sigma2 = sigma * sigma;
        
        // Partial derivatives
        T dSigma_dr = 2.0 * r;
        T dSigma_dtheta = -2.0 * a2 * sin_theta * cos_theta;
        T dDelta_dr = 2.0 * r - 2.0 * M_;
        
        // This is a simplified version - full expressions are more complex
        // Γ^t components
        Gamma[0][0][1] = M_ * (r2 - a2 * cos2_theta) / (sigma2 * delta) * (r2 + a2);
        Gamma[0][1][0] = Gamma[0][0][1];
        
        Gamma[0][0][2] = -2.0 * M_ * a2 * r * sin_theta * cos_theta / sigma2;
        Gamma[0][2][0] = Gamma[0][0][2];
        
        Gamma[0][1][3] = -a_ * M_ * sin2_theta * (r2 - a2 * cos2_theta) / (sigma2 * delta) * 
                         (r2 + a2 + 2.0 * M_ * r);
        Gamma[0][3][1] = Gamma[0][1][3];
        
        // Γ^r components
        Gamma[1][0][0] = M_ * delta * (r2 - a2 * cos2_theta) / (sigma2 * sigma);
        
        Gamma[1][1][1] = (r * delta - M_ * (r2 - a2)) / (sigma * delta);
        
        Gamma[1][2][2] = -r * delta / sigma;
        
        Gamma[1][3][3] = -delta * sin2_theta * (r - M_ * a2 * sin2_theta / sigma) / sigma;
        
        // Γ^θ components
        Gamma[2][1][2] = r / sigma;
        Gamma[2][2][1] = Gamma[2][1][2];
        
        Gamma[2][3][3] = -sin_theta * cos_theta * 
                         ((r2 + a2) * (r2 + a2) / sigma + a2 * sin2_theta) / sigma;
        
        // Γ^φ components
        Gamma[3][1][3] = (r * sigma - M_ * (r2 - a2 * cos2_theta)) / (sigma2 * delta);
        Gamma[3][3][1] = Gamma[3][1][3];
        
        Gamma[3][2][3] = cos_theta / sin_theta + 2.0 * a2 * sin_theta * cos_theta / sigma;
        Gamma[3][3][2] = Gamma[3][2][3];
        
        return Gamma;
    }
    
    /**
     * @brief Compute constants of motion for geodesic
     * @param position 4-position
     * @param velocity 4-velocity
     * @return {E, L, Q} energy, angular momentum, Carter constant
     */
    std::array<T, 3> constants_of_motion(const Vec4<T>& pos, const Vec4<T>& vel) const {
        auto metric = compute(pos.r, pos.theta);
        
        // Energy E = -p_t
        T E = 0;
        for (int nu = 0; nu < 4; ++nu) {
            E -= metric.g[0][nu] * vel[nu];
        }
        
        // Angular momentum L = p_φ
        T L = 0;
        for (int nu = 0; nu < 4; ++nu) {
            L += metric.g[3][nu] * vel[nu];
        }
        
        // Carter constant Q
        T sigma = Sigma(pos.r, pos.theta);
        T p_theta = metric.g[2][2] * vel.theta;
        T cos_theta = std::cos(pos.theta);
        T Q = p_theta * p_theta + cos_theta * cos_theta * 
              (a_ * a_ * (1.0 - E * E) + L * L / (std::sin(pos.theta) * std::sin(pos.theta)));
        
        return {E, L, Q};
    }
};

/**
 * @brief Kerr-Schild coordinates (horizon-penetrating)
 * 
 * These coordinates are regular at the event horizon, unlike Boyer-Lindquist.
 */
template<typename T = double>
class KerrSchildMetric : public KerrMetric<T> {
public:
    using KerrMetric<T>::KerrMetric;
    
    /**
     * @brief Transform from Boyer-Lindquist to Kerr-Schild coordinates
     */
    static Vec4<T> BL_to_KS(const Vec4<T>& bl, T M, T a) {
        Vec4<T> ks;
        T r = bl.r;
        T delta = r * r - 2.0 * M * r + a * a;
        
        ks.t = bl.t + 2.0 * M * r / delta;
        ks.r = bl.r;
        ks.theta = bl.theta;
        ks.phi = bl.phi + a / delta;
        
        return ks;
    }
    
    /**
     * @brief Transform from Kerr-Schild to Boyer-Lindquist coordinates
     */
    static Vec4<T> KS_to_BL(const Vec4<T>& ks, T M, T a) {
        Vec4<T> bl;
        T r = ks.r;
        T delta = r * r - 2.0 * M * r + a * a;
        
        bl.t = ks.t - 2.0 * M * r / delta;
        bl.r = ks.r;
        bl.theta = ks.theta;
        bl.phi = ks.phi - a / delta;
        
        return bl;
    }
};

// Factory function to create metrics
template<typename T = double>
std::unique_ptr<Metric<T>> create_metric(const std::string& type, T mass, T spin = 0.0) {
    if (type == "schwarzschild") {
        return std::make_unique<SchwarzschildMetric<T>>(mass);
    } else if (type == "kerr" || type == "kerr-bl") {
        return std::make_unique<KerrMetric<T>>(mass, spin);
    } else if (type == "kerr-schild" || type == "kerr-ks") {
        return std::make_unique<KerrSchildMetric<T>>(mass, spin);
    } else {
        throw std::invalid_argument("Unknown metric type: " + type);
    }
}

} // namespace blackhole

#endif // BLACKHOLE_METRIC_HPP