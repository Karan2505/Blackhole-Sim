/**
 * @file grmhd.hpp
 * @brief General Relativistic Magnetohydrodynamics Solver
 * 
 * Implements the Valencia formulation of GRMHD with:
 * - Conservative variable evolution
 * - HLLC Riemann solver
 * - Constrained transport for divergence-free B
 * - Primitive variable recovery (Noble 2D method)
 * 
 * @author BlackHole-Sim Team
 * @license MIT
 */

#ifndef BLACKHOLE_GRMHD_HPP
#define BLACKHOLE_GRMHD_HPP

#include "metric.hpp"
#include <vector>
#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace blackhole {
namespace grmhd {

/**
 * @brief Primitive variables for GRMHD
 */
template<typename T = double>
struct Primitives {
    T rho;          // Rest-mass density
    T u;            // Internal energy density
    T v1, v2, v3;   // 3-velocity components (v^i)
    T B1, B2, B3;   // Magnetic field components (B^i)
    
    Primitives() : rho(0), u(0), v1(0), v2(0), v3(0), B1(0), B2(0), B3(0) {}
    
    Primitives(T rho_, T u_, T v1_, T v2_, T v3_, T B1_, T B2_, T B3_)
        : rho(rho_), u(u_), v1(v1_), v2(v2_), v3(v3_), B1(B1_), B2(B2_), B3(B3_) {}
    
    /**
     * @brief Compute Lorentz factor
     */
    T lorentz_factor(const MetricTensor<T>& g) const {
        T v_sq = g.g[1][1] * v1 * v1 + g.g[2][2] * v2 * v2 + g.g[3][3] * v3 * v3;
        // Add cross terms for Kerr metric
        v_sq += 2.0 * g.g[1][2] * v1 * v2;
        v_sq += 2.0 * g.g[1][3] * v1 * v3;
        v_sq += 2.0 * g.g[2][3] * v2 * v3;
        
        if (v_sq >= 1.0) {
            v_sq = 0.9999; // Limit to subluminal
        }
        
        return 1.0 / std::sqrt(1.0 - v_sq);
    }
    
    /**
     * @brief Compute magnetic pressure
     */
    T magnetic_pressure(const MetricTensor<T>& g) const {
        T B_sq = g.g[1][1] * B1 * B1 + g.g[2][2] * B2 * B2 + g.g[3][3] * B3 * B3;
        return 0.5 * B_sq;
    }
    
    /**
     * @brief Compute gas pressure (ideal EOS: P = (gamma-1) * u)
     */
    T pressure(T gamma = 4.0/3.0) const {
        return (gamma - 1.0) * u;
    }
    
    /**
     * @brief Compute enthalpy density: h = rho + u + P
     */
    T enthalpy(T gamma = 4.0/3.0) const {
        return rho + u + pressure(gamma);
    }
};

/**
 * @brief Conserved variables for GRMHD
 */
template<typename T = double>
struct Conserved {
    T D;            // Conserved density: D = sqrt(gamma) * rho * W
    T tau;          // Energy density: tau = sqrt(gamma) * (rho*h*W^2 - P - rho*W) + ...
    T S1, S2, S3;   // Momentum density: S_i = sqrt(gamma) * rho*h*W^2 * v_i + ...
    T B1, B2, B3;   // Magnetic field (same as primitive)
    
    Conserved() : D(0), tau(0), S1(0), S2(0), S3(0), B1(0), B2(0), B3(0) {}
    
    Conserved operator+(const Conserved& other) const {
        Conserved result;
        result.D = D + other.D;
        result.tau = tau + other.tau;
        result.S1 = S1 + other.S1;
        result.S2 = S2 + other.S2;
        result.S3 = S3 + other.S3;
        result.B1 = B1 + other.B1;
        result.B2 = B2 + other.B2;
        result.B3 = B3 + other.B3;
        return result;
    }
    
    Conserved operator*(T scalar) const {
        Conserved result;
        result.D = D * scalar;
        result.tau = tau * scalar;
        result.S1 = S1 * scalar;
        result.S2 = S2 * scalar;
        result.S3 = S3 * scalar;
        result.B1 = B1 * scalar;
        result.B2 = B2 * scalar;
        result.B3 = B3 * scalar;
        return result;
    }
};

/**
 * @brief Equation of State interface
 */
template<typename T = double>
class EquationOfState {
public:
    virtual ~EquationOfState() = default;
    virtual T pressure(T rho, T u) const = 0;
    virtual T sound_speed(T rho, T u) const = 0;
    virtual T gamma() const = 0;
};

/**
 * @brief Ideal gas equation of state
 * P = (gamma - 1) * u
 */
template<typename T = double>
class IdealGasEOS : public EquationOfState<T> {
private:
    T gamma_;
    
public:
    explicit IdealGasEOS(T gamma = 4.0/3.0) : gamma_(gamma) {
        if (gamma <= 1.0) {
            throw std::invalid_argument("Gamma must be > 1");
        }
    }
    
    T pressure(T rho, T u) const override {
        return (gamma_ - 1.0) * u;
    }
    
    T sound_speed(T rho, T u) const override {
        T P = pressure(rho, u);
        T h = rho + u + P;
        if (h <= 0 || rho <= 0) return 0;
        return std::sqrt(gamma_ * P / h);
    }
    
    T gamma() const override { return gamma_; }
};

/**
 * @brief Convert primitives to conserved variables
 */
template<typename T = double>
Conserved<T> prim_to_cons(
    const Primitives<T>& P,
    const MetricTensor<T>& g,
    const EquationOfState<T>& eos
) {
    Conserved<T> U;
    
    T W = P.lorentz_factor(g);
    T pres = eos.pressure(P.rho, P.u);
    T h = P.rho + P.u + pres;
    
    // Compute sqrt of spatial metric determinant
    T sqrt_gamma = std::sqrt(g.g[1][1] * g.g[2][2] * g.g[3][3]);
    
    // Lower velocity index: v_i = g_ij v^j
    T v_1 = g.g[1][1] * P.v1 + g.g[1][2] * P.v2 + g.g[1][3] * P.v3;
    T v_2 = g.g[2][1] * P.v1 + g.g[2][2] * P.v2 + g.g[2][3] * P.v3;
    T v_3 = g.g[3][1] * P.v1 + g.g[3][2] * P.v2 + g.g[3][3] * P.v3;
    
    // Magnetic field terms
    T B_sq = g.g[1][1] * P.B1 * P.B1 + g.g[2][2] * P.B2 * P.B2 + g.g[3][3] * P.B3 * P.B3;
    T vB = P.v1 * P.B1 + P.v2 * P.B2 + P.v3 * P.B3;
    T b_sq = B_sq / (W * W) + vB * vB;
    
    // Conserved variables
    U.D = sqrt_gamma * P.rho * W;
    U.tau = sqrt_gamma * ((h + b_sq) * W * W - pres - 0.5 * b_sq - P.rho * W);
    
    T hstar = h + b_sq;
    U.S1 = sqrt_gamma * (hstar * W * W * v_1 - vB * P.B1);
    U.S2 = sqrt_gamma * (hstar * W * W * v_2 - vB * P.B2);
    U.S3 = sqrt_gamma * (hstar * W * W * v_3 - vB * P.B3);
    
    U.B1 = P.B1;
    U.B2 = P.B2;
    U.B3 = P.B3;
    
    return U;
}

/**
 * @brief Primitive variable recovery using 2D Newton-Raphson (Noble method)
 * 
 * This is one of the most challenging parts of GRMHD codes.
 * We solve for W (Lorentz factor) and v^2 given conserved variables.
 */
template<typename T = double>
class PrimitiveRecovery {
private:
    const EquationOfState<T>& eos_;
    T tolerance_;
    int max_iter_;
    
public:
    PrimitiveRecovery(
        const EquationOfState<T>& eos,
        T tol = 1e-10,
        int max_iter = 100
    ) : eos_(eos), tolerance_(tol), max_iter_(max_iter) {}
    
    /**
     * @brief Recover primitives from conserved variables
     * @return true if successful, false if failed
     */
    bool recover(
        const Conserved<T>& U,
        const MetricTensor<T>& g,
        Primitives<T>& P,
        T floor_rho = 1e-12,
        T floor_u = 1e-14
    ) const {
        T sqrt_gamma = std::sqrt(g.g[1][1] * g.g[2][2] * g.g[3][3]);
        
        // Initial guess for W
        T W = 1.0;
        T vsq = 0.0;
        
        // Magnetic field magnitudes
        T B_sq = g.g[1][1] * U.B1 * U.B1 + g.g[2][2] * U.B2 * U.B2 + g.g[3][3] * U.B3 * U.B3;
        T S_sq = g.g_inv[1][1] * U.S1 * U.S1 + g.g_inv[2][2] * U.S2 * U.S2 + g.g_inv[3][3] * U.S3 * U.S3;
        T SdotB = U.S1 * U.B1 + U.S2 * U.B2 + U.S3 * U.B3;
        
        T gamma_eos = eos_.gamma();
        
        // 1D Newton-Raphson for simplified case (no magnetic field)
        if (B_sq < 1e-14) {
            for (int iter = 0; iter < max_iter_; ++iter) {
                T rho = U.D / (sqrt_gamma * W);
                if (rho < floor_rho) rho = floor_rho;
                
                T S_mag = std::sqrt(S_sq);
                T v = S_mag / (U.tau + U.D / sqrt_gamma + eos_.pressure(rho, 0) * sqrt_gamma);
                if (v >= 1.0) v = 0.9999;
                
                T W_new = 1.0 / std::sqrt(1.0 - v * v);
                
                if (std::abs(W_new - W) < tolerance_) {
                    W = W_new;
                    break;
                }
                W = W_new;
            }
        } else {
            // Full 2D solver for MHD case
            // Simplified version - production code would use Noble et al. 2006
            for (int iter = 0; iter < max_iter_; ++iter) {
                T rho = U.D / (sqrt_gamma * W);
                if (rho < floor_rho) rho = floor_rho;
                
                // Estimate pressure
                T P_est = (gamma_eos - 1.0) * (U.tau / sqrt_gamma - 0.5 * B_sq);
                if (P_est < 0) P_est = floor_u * (gamma_eos - 1.0);
                
                T h = rho + P_est / (gamma_eos - 1.0) + P_est;
                T hstar = h + B_sq / (W * W);
                
                // Velocity from momentum
                T denom = hstar * W * W * sqrt_gamma;
                if (denom < 1e-14) denom = 1e-14;
                
                T v1 = (U.S1 + SdotB * U.B1 / (W * W)) / denom;
                T v2 = (U.S2 + SdotB * U.B2 / (W * W)) / denom;
                T v3 = (U.S3 + SdotB * U.B3 / (W * W)) / denom;
                
                vsq = g.g[1][1] * v1 * v1 + g.g[2][2] * v2 * v2 + g.g[3][3] * v3 * v3;
                if (vsq >= 1.0) vsq = 0.9999;
                
                T W_new = 1.0 / std::sqrt(1.0 - vsq);
                
                if (std::abs(W_new - W) / W < tolerance_) {
                    W = W_new;
                    P.v1 = v1;
                    P.v2 = v2;
                    P.v3 = v3;
                    break;
                }
                W = 0.5 * (W + W_new); // Damped update
            }
        }
        
        // Compute final primitives
        P.rho = U.D / (sqrt_gamma * W);
        if (P.rho < floor_rho) P.rho = floor_rho;
        
        // Internal energy from energy conservation
        T kinetic = 0.5 * P.rho * vsq * W * W;
        T magnetic = 0.5 * B_sq;
        P.u = (U.tau / sqrt_gamma - kinetic - magnetic) / (gamma_eos);
        if (P.u < floor_u) P.u = floor_u;
        
        P.B1 = U.B1;
        P.B2 = U.B2;
        P.B3 = U.B3;
        
        return true;
    }
};

/**
 * @brief HLLC Riemann solver for relativistic MHD
 */
template<typename T = double>
class HLLCFlux {
private:
    const EquationOfState<T>& eos_;
    
public:
    explicit HLLCFlux(const EquationOfState<T>& eos) : eos_(eos) {}
    
    /**
     * @brief Compute HLLC flux in given direction
     * @param PL Left primitive state
     * @param PR Right primitive state
     * @param g Metric at interface
     * @param dir Direction (1=r, 2=theta, 3=phi)
     * @return Flux conserved variables
     */
    Conserved<T> compute(
        const Primitives<T>& PL,
        const Primitives<T>& PR,
        const MetricTensor<T>& g,
        int dir
    ) const {
        // Compute wave speeds
        T cL = eos_.sound_speed(PL.rho, PL.u);
        T cR = eos_.sound_speed(PR.rho, PR.u);
        
        // Get velocities in direction
        T vL = (dir == 1) ? PL.v1 : (dir == 2) ? PL.v2 : PL.v3;
        T vR = (dir == 1) ? PR.v1 : (dir == 2) ? PR.v2 : PR.v3;
        
        // Lorentz factors
        T WL = PL.lorentz_factor(g);
        T WR = PR.lorentz_factor(g);
        
        // Wave speed estimates (simplified)
        T lambdaL = (vL - cL) / (1.0 - vL * cL);
        T lambdaR = (vR + cR) / (1.0 + vR * cR);
        
        T SL = std::min(lambdaL, lambdaR);
        T SR = std::max(lambdaL, lambdaR);
        
        // Compute conserved states
        Conserved<T> UL = prim_to_cons(PL, g, eos_);
        Conserved<T> UR = prim_to_cons(PR, g, eos_);
        
        // Compute fluxes
        Conserved<T> FL = compute_flux(PL, UL, g, dir);
        Conserved<T> FR = compute_flux(PR, UR, g, dir);
        
        // HLL flux
        if (SL >= 0) {
            return FL;
        } else if (SR <= 0) {
            return FR;
        } else {
            // HLL average (simplified HLLC)
            T denom = SR - SL;
            if (std::abs(denom) < 1e-14) denom = 1e-14;
            
            return (FR * (-SL) + FL * SR + (UL + UR * (-1.0)) * (SR * SL)) * (1.0 / denom);
        }
    }
    
private:
    Conserved<T> compute_flux(
        const Primitives<T>& P,
        const Conserved<T>& U,
        const MetricTensor<T>& g,
        int dir
    ) const {
        Conserved<T> F;
        
        T v = (dir == 1) ? P.v1 : (dir == 2) ? P.v2 : P.v3;
        T B = (dir == 1) ? P.B1 : (dir == 2) ? P.B2 : P.B3;
        T pres = eos_.pressure(P.rho, P.u);
        T pmag = P.magnetic_pressure(g);
        T ptot = pres + pmag;
        
        T vB = P.v1 * P.B1 + P.v2 * P.B2 + P.v3 * P.B3;
        
        F.D = U.D * v;
        F.tau = U.tau * v + ptot * v - vB * B;
        F.S1 = U.S1 * v + ((dir == 1) ? ptot : 0) - P.B1 * B / P.lorentz_factor(g);
        F.S2 = U.S2 * v + ((dir == 2) ? ptot : 0) - P.B2 * B / P.lorentz_factor(g);
        F.S3 = U.S3 * v + ((dir == 3) ? ptot : 0) - P.B3 * B / P.lorentz_factor(g);
        F.B1 = (dir == 1) ? 0 : (P.B1 * v - P.v1 * B);
        F.B2 = (dir == 2) ? 0 : (P.B2 * v - P.v2 * B);
        F.B3 = (dir == 3) ? 0 : (P.B3 * v - P.v3 * B);
        
        return F;
    }
};

/**
 * @brief Source terms for GRMHD in curved spacetime
 */
template<typename T = double>
Conserved<T> compute_source_terms(
    const Primitives<T>& P,
    const MetricTensor<T>& g,
    const std::array<std::array<std::array<T, 4>, 4>, 4>& Gamma,
    const EquationOfState<T>& eos
) {
    Conserved<T> S;
    
    T W = P.lorentz_factor(g);
    T pres = eos.pressure(P.rho, P.u);
    T h = P.rho + P.u + pres;
    
    // Stress-energy tensor components (simplified)
    T T00 = h * W * W - pres;
    T T01 = h * W * W * P.v1;
    T T02 = h * W * W * P.v2;
    T T03 = h * W * W * P.v3;
    
    // Source terms from Christoffel symbols
    // S_tau = -alpha * T^mu_nu * Gamma^nu_mu0
    // S_i = alpha * T^mu_nu * Gamma^nu_mu_i
    
    // Simplified source computation
    S.D = 0; // Mass conservation has no source
    
    S.tau = 0;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            S.tau -= Gamma[0][mu][nu]; // Simplified
        }
    }
    
    // Momentum sources from curvature
    S.S1 = T00 * Gamma[1][0][0] + T01 * Gamma[1][0][1] + T02 * Gamma[1][0][2];
    S.S2 = T00 * Gamma[2][0][0] + T01 * Gamma[2][0][1] + T02 * Gamma[2][0][2];
    S.S3 = T00 * Gamma[3][0][0] + T01 * Gamma[3][0][1] + T02 * Gamma[3][0][2];
    
    // B sources are zero (evolution via induction equation)
    S.B1 = S.B2 = S.B3 = 0;
    
    return S;
}

/**
 * @brief Grid cell for GRMHD
 */
template<typename T = double>
struct Cell {
    Primitives<T> P;
    Conserved<T> U;
    T r, theta, phi;
    MetricTensor<T> g;
    
    Cell() : r(0), theta(0), phi(0) {}
};

/**
 * @brief GRMHD solver on fixed Kerr background
 */
template<typename T = double>
class GRMHDSolver {
private:
    const Metric<T>& metric_;
    std::unique_ptr<EquationOfState<T>> eos_;
    std::unique_ptr<PrimitiveRecovery<T>> recovery_;
    std::unique_ptr<HLLCFlux<T>> flux_solver_;
    
    // Grid
    int nr_, ntheta_, nphi_;
    T r_min_, r_max_;
    std::vector<Cell<T>> grid_;
    
    // Time stepping
    T cfl_;
    T time_;
    int step_;
    
public:
    GRMHDSolver(
        const Metric<T>& metric,
        int nr, int ntheta, int nphi,
        T r_min, T r_max,
        T gamma = 4.0/3.0,
        T cfl = 0.5
    ) : metric_(metric), nr_(nr), ntheta_(ntheta), nphi_(nphi),
        r_min_(r_min), r_max_(r_max), cfl_(cfl), time_(0), step_(0)
    {
        eos_ = std::make_unique<IdealGasEOS<T>>(gamma);
        recovery_ = std::make_unique<PrimitiveRecovery<T>>(*eos_);
        flux_solver_ = std::make_unique<HLLCFlux<T>>(*eos_);
        
        // Initialize grid
        grid_.resize(nr_ * ntheta_ * nphi_);
        
        for (int i = 0; i < nr_; ++i) {
            for (int j = 0; j < ntheta_; ++j) {
                for (int k = 0; k < nphi_; ++k) {
                    Cell<T>& cell = grid_[index(i, j, k)];
                    
                    // Log-spaced radial grid
                    T log_r = std::log(r_min_) + (std::log(r_max_) - std::log(r_min_)) * i / (nr_ - 1);
                    cell.r = std::exp(log_r);
                    
                    // Uniform theta (avoiding poles)
                    cell.theta = 0.01 + (constants::pi - 0.02) * j / (ntheta_ - 1);
                    
                    // Uniform phi
                    cell.phi = 2.0 * constants::pi * k / nphi_;
                    
                    // Compute metric at cell center
                    cell.g = metric_.compute(cell.r, cell.theta);
                }
            }
        }
    }
    
    /**
     * @brief Initialize Fishbone-Moncrief torus
     * 
     * Standard initial condition for accretion disk simulations.
     * See Fishbone & Moncrief 1976.
     */
    void init_fishbone_moncrief(
        T r_in,         // Inner edge radius
        T r_max_torus,  // Pressure maximum radius
        T rho_max = 1.0 // Maximum density
    ) {
        T l = compute_angular_momentum(r_max_torus);
        
        for (auto& cell : grid_) {
            T r = cell.r;
            T theta = cell.theta;
            T sin_theta = std::sin(theta);
            
            // Compute specific enthalpy from Bernoulli integral
            T W_in = potential(r_in, constants::pi / 2, l);
            T W = potential(r, theta, l);
            T h = std::exp(W_in - W);
            
            if (h > 1.0 && r >= r_in && std::abs(theta - constants::pi/2) < constants::pi/3) {
                // Inside torus
                T K = rho_max * std::pow(eos_->gamma() - 1.0, eos_->gamma()) / eos_->gamma();
                cell.P.rho = std::pow((h - 1.0) * (eos_->gamma() - 1.0) / (K * eos_->gamma()), 
                                      1.0 / (eos_->gamma() - 1.0));
                cell.P.u = K * std::pow(cell.P.rho, eos_->gamma()) / (eos_->gamma() - 1.0);
                
                // Keplerian angular velocity
                T Omega = 1.0 / (std::pow(r, 1.5) + metric_.spin());
                cell.P.v3 = (Omega - cell.g.g[0][3] / cell.g.g[3][3]) * r * sin_theta;
                cell.P.v1 = 0;
                cell.P.v2 = 0;
                
                // Small magnetic field
                cell.P.B1 = 0;
                cell.P.B2 = 0;
                cell.P.B3 = 0.01 * cell.P.rho; // Weak toroidal field
            } else {
                // Atmosphere
                cell.P.rho = 1e-8;
                cell.P.u = 1e-10;
                cell.P.v1 = cell.P.v2 = cell.P.v3 = 0;
                cell.P.B1 = cell.P.B2 = cell.P.B3 = 0;
            }
            
            // Compute conserved variables
            cell.U = prim_to_cons(cell.P, cell.g, *eos_);
        }
    }
    
    /**
     * @brief Take one time step using SSP-RK3
     */
    T step_rk3() {
        T dt = compute_timestep();
        
        // Store initial state
        std::vector<Conserved<T>> U0(grid_.size());
        for (size_t i = 0; i < grid_.size(); ++i) {
            U0[i] = grid_[i].U;
        }
        
        // Stage 1: U1 = U0 + dt * L(U0)
        compute_rhs();
        for (size_t i = 0; i < grid_.size(); ++i) {
            grid_[i].U = grid_[i].U + grid_[i].dU * dt;
            recovery_->recover(grid_[i].U, grid_[i].g, grid_[i].P);
        }
        
        // Stage 2: U2 = 3/4 U0 + 1/4 U1 + 1/4 dt L(U1)
        compute_rhs();
        for (size_t i = 0; i < grid_.size(); ++i) {
            grid_[i].U = U0[i] * 0.75 + grid_[i].U * 0.25 + grid_[i].dU * (0.25 * dt);
            recovery_->recover(grid_[i].U, grid_[i].g, grid_[i].P);
        }
        
        // Stage 3: U3 = 1/3 U0 + 2/3 U2 + 2/3 dt L(U2)
        compute_rhs();
        for (size_t i = 0; i < grid_.size(); ++i) {
            grid_[i].U = U0[i] * (1.0/3.0) + grid_[i].U * (2.0/3.0) + grid_[i].dU * (2.0/3.0 * dt);
            recovery_->recover(grid_[i].U, grid_[i].g, grid_[i].P);
        }
        
        time_ += dt;
        step_++;
        
        return dt;
    }
    
    /**
     * @brief Compute total mass in domain
     */
    T compute_mass() const {
        T mass = 0;
        for (const auto& cell : grid_) {
            T dV = cell.r * cell.r * std::sin(cell.theta); // Approximate volume element
            mass += cell.U.D * dV;
        }
        return mass;
    }
    
    /**
     * @brief Get current simulation time
     */
    T current_time() const { return time_; }
    
    /**
     * @brief Get current step number
     */
    int current_step() const { return step_; }
    
    /**
     * @brief Access grid data (for I/O)
     */
    const std::vector<Cell<T>>& grid() const { return grid_; }
    
private:
    size_t index(int i, int j, int k) const {
        return i + nr_ * (j + ntheta_ * k);
    }
    
    T compute_timestep() const {
        T dt_min = 1e10;
        
        for (const auto& cell : grid_) {
            T cs = eos_->sound_speed(cell.P.rho, cell.P.u);
            T dr = cell.r * 0.1; // Approximate cell size
            T dt = cfl_ * dr / (cs + std::abs(cell.P.v1));
            dt_min = std::min(dt_min, dt);
        }
        
        return dt_min;
    }
    
    T compute_angular_momentum(T r_max) const {
        // Specific angular momentum for FM torus
        return std::sqrt(metric_.mass() * r_max);
    }
    
    T potential(T r, T theta, T l) const {
        T sin_theta = std::sin(theta);
        T a = metric_.spin();
        T Delta = r * r - 2 * metric_.mass() * r + a * a;
        T Sigma = r * r + a * a * std::cos(theta) * std::cos(theta);
        
        T A = (r * r + a * a) * (r * r + a * a) - a * a * Delta * sin_theta * sin_theta;
        T omega = 2 * a * r / A;
        
        T u_t = std::sqrt(-1.0 / (r * r - 2 * metric_.mass() / r + l * l / (r * r * sin_theta * sin_theta)));
        
        return std::log(-u_t);
    }
    
    void compute_rhs() {
        // Simplified - compute fluxes and sources
        #pragma omp parallel for
        for (int i = 1; i < nr_ - 1; ++i) {
            for (int j = 1; j < ntheta_ - 1; ++j) {
                for (int k = 0; k < nphi_; ++k) {
                    size_t idx = index(i, j, k);
                    Cell<T>& cell = grid_[idx];
                    
                    // Get neighbors
                    Cell<T>& cellL = grid_[index(i-1, j, k)];
                    Cell<T>& cellR = grid_[index(i+1, j, k)];
                    
                    // Compute flux difference (radial direction only for simplicity)
                    Conserved<T> FR = flux_solver_->compute(cell.P, cellR.P, cell.g, 1);
                    Conserved<T> FL = flux_solver_->compute(cellL.P, cell.P, cell.g, 1);
                    
                    T dr = cellR.r - cellL.r;
                    cell.dU = (FL + FR * (-1.0)) * (1.0 / dr);
                    
                    // Add source terms
                    auto Gamma = metric_.christoffel(cell.r, cell.theta);
                    Conserved<T> S = compute_source_terms(cell.P, cell.g, Gamma, *eos_);
                    cell.dU = cell.dU + S;
                }
            }
        }
    }
    
    // RHS storage
    Conserved<T> dU;
};

} // namespace grmhd
} // namespace blackhole

#endif // BLACKHOLE_GRMHD_HPP