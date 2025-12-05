/**
 * @file test_metric.cpp
 * @brief Unit tests for spacetime metric implementations
 * 
 * Tests include:
 * - Schwarzschild metric properties
 * - Kerr metric properties
 * - ISCO radius validation
 * - Photon sphere validation
 * - Christoffel symbol symmetry
 * 
 * @author BlackHole-Sim Team
 * @license MIT
 */

#include <gtest/gtest.h>
#include <cmath>
#include "blackhole/metric.hpp"

using namespace blackhole;

// ============================================================================
// Schwarzschild Metric Tests
// ============================================================================

class SchwarzschildMetricTest : public ::testing::Test {
protected:
    SchwarzschildMetric<double> metric{1.0};  // Unit mass
};

TEST_F(SchwarzschildMetricTest, HorizonRadius) {
    // r_s = 2M for Schwarzschild
    EXPECT_DOUBLE_EQ(metric.horizon_radius(), 2.0);
}

TEST_F(SchwarzschildMetricTest, ISCORadius) {
    // r_ISCO = 6M for Schwarzschild
    EXPECT_DOUBLE_EQ(metric.isco_radius(), 6.0);
}

TEST_F(SchwarzschildMetricTest, PhotonSphereRadius) {
    // r_photon = 3M for Schwarzschild
    EXPECT_DOUBLE_EQ(metric.photon_sphere_radius(), 3.0);
}

TEST_F(SchwarzschildMetricTest, MetricSignature) {
    // Test metric signature (-,+,+,+) at r = 10M
    auto g = metric.compute(10.0, M_PI / 2.0);
    
    EXPECT_LT(g.g[0][0], 0);  // g_tt < 0
    EXPECT_GT(g.g[1][1], 0);  // g_rr > 0
    EXPECT_GT(g.g[2][2], 0);  // g_θθ > 0
    EXPECT_GT(g.g[3][3], 0);  // g_φφ > 0
}

TEST_F(SchwarzschildMetricTest, MetricDeterminant) {
    // |g| = -r^4 sin^2(θ) for Schwarzschild
    double r = 5.0;
    double theta = M_PI / 3.0;
    auto g = metric.compute(r, theta);
    
    double expected_det = -std::pow(r, 4) * std::pow(std::sin(theta), 2);
    EXPECT_NEAR(g.determinant(), expected_det, 1e-10);
}

TEST_F(SchwarzschildMetricTest, AsymptoticFlatness) {
    // At large r, metric approaches Minkowski
    double r = 1000.0;
    auto g = metric.compute(r, M_PI / 2.0);
    
    EXPECT_NEAR(g.g[0][0], -1.0, 0.01);
    EXPECT_NEAR(g.g[1][1], 1.0, 0.01);
}

TEST_F(SchwarzschildMetricTest, ChristoffelSymmetry) {
    // Γ^α_μν = Γ^α_νμ (symmetric in lower indices)
    double r = 5.0;
    double theta = M_PI / 3.0;
    auto Gamma = metric.christoffel(r, theta);
    
    for (int alpha = 0; alpha < 4; alpha++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu; nu < 4; nu++) {
                EXPECT_NEAR(Gamma[alpha][mu][nu], Gamma[alpha][nu][mu], 1e-14)
                    << "Asymmetry at [" << alpha << "][" << mu << "][" << nu << "]";
            }
        }
    }
}

TEST_F(SchwarzschildMetricTest, InverseMetric) {
    // g * g^(-1) = identity
    double r = 7.0;
    double theta = M_PI / 4.0;
    auto g = metric.compute(r, theta);
    
    // Check diagonal elements
    EXPECT_NEAR(g.g[0][0] * g.g_inv[0][0], 1.0, 1e-12);
    EXPECT_NEAR(g.g[1][1] * g.g_inv[1][1], 1.0, 1e-12);
    EXPECT_NEAR(g.g[2][2] * g.g_inv[2][2], 1.0, 1e-12);
    EXPECT_NEAR(g.g[3][3] * g.g_inv[3][3], 1.0, 1e-12);
}

// ============================================================================
// Kerr Metric Tests
// ============================================================================

class KerrMetricTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Near-extremal Kerr (like Interstellar's Gargantua)
        metric_extreme = std::make_unique<KerrMetric<double>>(1.0, 0.998);
        // Moderate spin
        metric_moderate = std::make_unique<KerrMetric<double>>(1.0, 0.5);
        // Slow spin
        metric_slow = std::make_unique<KerrMetric<double>>(1.0, 0.1);
    }
    
    std::unique_ptr<KerrMetric<double>> metric_extreme;
    std::unique_ptr<KerrMetric<double>> metric_moderate;
    std::unique_ptr<KerrMetric<double>> metric_slow;
};

TEST_F(KerrMetricTest, HorizonRadiusExtreme) {
    // r_+ = M + sqrt(M^2 - a^2) for a = 0.998M
    double expected = 1.0 + std::sqrt(1.0 - 0.998 * 0.998);
    EXPECT_NEAR(metric_extreme->horizon_radius(), expected, 1e-10);
}

TEST_F(KerrMetricTest, HorizonRadiusSchwarzschild) {
    // As a → 0, r_+ → 2M (Schwarzschild limit)
    KerrMetric<double> nearly_schwarzschild(1.0, 0.001);
    EXPECT_NEAR(nearly_schwarzschild.horizon_radius(), 2.0, 0.01);
}

TEST_F(KerrMetricTest, ISCODecreaseWithSpin) {
    // ISCO moves inward for prograde orbits as spin increases
    double isco_slow = metric_slow->isco_radius();
    double isco_mod = metric_moderate->isco_radius();
    double isco_ext = metric_extreme->isco_radius();
    
    EXPECT_GT(isco_slow, isco_mod);
    EXPECT_GT(isco_mod, isco_ext);
}

TEST_F(KerrMetricTest, ISCOLimit) {
    // For a → M, r_ISCO → M (prograde)
    EXPECT_LT(metric_extreme->isco_radius(), 1.5);  // Should be close to M
}

TEST_F(KerrMetricTest, PhotonSphereExists) {
    // Photon sphere should exist for all Kerr black holes
    EXPECT_GT(metric_extreme->photon_sphere_radius(), metric_extreme->horizon_radius());
    EXPECT_GT(metric_moderate->photon_sphere_radius(), metric_moderate->horizon_radius());
}

TEST_F(KerrMetricTest, ErgosphereAtEquator) {
    // At equator (θ = π/2), ergosphere is at r = 2M (same as Schwarzschild horizon)
    double r_ergo = metric_moderate->ergosphere_radius(M_PI / 2.0);
    EXPECT_DOUBLE_EQ(r_ergo, 2.0);
}

TEST_F(KerrMetricTest, ErgosphereAtPoles) {
    // At poles (θ = 0 or π), ergosphere touches horizon
    double r_ergo_pole = metric_moderate->ergosphere_radius(0.01);  // Near pole
    double r_horizon = metric_moderate->horizon_radius();
    EXPECT_NEAR(r_ergo_pole, r_horizon, 0.1);
}

TEST_F(KerrMetricTest, FrameDragging) {
    // g_tφ ≠ 0 for rotating black hole (frame-dragging)
    auto g = metric_moderate->compute(5.0, M_PI / 2.0);
    EXPECT_NE(g.g[0][3], 0.0);
}

TEST_F(KerrMetricTest, NoFrameDraggingOnAxis) {
    // g_tφ → 0 as θ → 0 or π (on rotation axis)
    auto g_axis = metric_moderate->compute(5.0, 0.01);
    EXPECT_NEAR(g_axis.g[0][3], 0.0, 0.01);
}

TEST_F(KerrMetricTest, HorizonAngularVelocity) {
    // Ω_H = a / (2 M r_+)
    double a = 0.5;
    double M = 1.0;
    double r_plus = M + std::sqrt(M * M - a * a);
    double expected_omega = a / (2.0 * M * r_plus);
    
    EXPECT_NEAR(metric_moderate->horizon_angular_velocity(), expected_omega, 1e-10);
}

TEST_F(KerrMetricTest, RetrogadeISCOLarger) {
    // Retrograde ISCO > Prograde ISCO
    double prograde = metric_moderate->isco_radius();
    double retrograde = metric_moderate->isco_radius_retrograde();
    
    EXPECT_GT(retrograde, prograde);
}

TEST_F(KerrMetricTest, ConstantsOfMotion) {
    // Test that constants are computed for a simple equatorial orbit
    Vec4<double> pos(0, 10.0, M_PI / 2.0, 0);
    
    // Circular orbit velocity
    double r = 10.0;
    double Omega = 1.0 / (std::pow(r, 1.5) + 0.5);  // Approximate
    Vec4<double> vel(1.0, 0, 0, Omega);
    
    auto [E, L, Q] = metric_moderate->constants_of_motion(pos, vel);
    
    // Energy should be positive for bound orbit
    EXPECT_GT(E, 0);
    // Angular momentum should be positive for prograde orbit
    EXPECT_GT(L, 0);
    // Carter constant should be zero for equatorial orbit
    EXPECT_NEAR(Q, 0, 1e-10);
}

// ============================================================================
// Metric Factory Tests
// ============================================================================

TEST(MetricFactoryTest, CreateSchwarzschild) {
    auto metric = create_metric<double>("schwarzschild", 2.0);
    EXPECT_DOUBLE_EQ(metric->mass(), 2.0);
    EXPECT_DOUBLE_EQ(metric->spin(), 0.0);
}

TEST(MetricFactoryTest, CreateKerr) {
    auto metric = create_metric<double>("kerr", 1.5, 0.7);
    EXPECT_DOUBLE_EQ(metric->mass(), 1.5);
    EXPECT_NEAR(metric->spin(), 0.7, 1e-10);
}

TEST(MetricFactoryTest, InvalidType) {
    EXPECT_THROW(create_metric<double>("invalid", 1.0), std::invalid_argument);
}

TEST(MetricFactoryTest, InvalidMass) {
    EXPECT_THROW(create_metric<double>("schwarzschild", -1.0), std::invalid_argument);
    EXPECT_THROW(create_metric<double>("schwarzschild", 0.0), std::invalid_argument);
}

TEST(MetricFactoryTest, InvalidSpin) {
    // Spin must be |a| < M for black hole
    EXPECT_THROW(create_metric<double>("kerr", 1.0, 1.0), std::invalid_argument);
    EXPECT_THROW(create_metric<double>("kerr", 1.0, 1.5), std::invalid_argument);
}

// ============================================================================
// Coordinate Transformation Tests
// ============================================================================

TEST(KerrSchildTest, BLtoKSAndBack) {
    double M = 1.0;
    double a = 0.5;
    
    Vec4<double> bl_original(10.0, 8.0, M_PI / 3.0, 1.5);
    
    // BL → KS → BL should recover original
    auto ks = KerrSchildMetric<double>::BL_to_KS(bl_original, M, a);
    auto bl_recovered = KerrSchildMetric<double>::KS_to_BL(ks, M, a);
    
    EXPECT_NEAR(bl_recovered.t, bl_original.t, 1e-10);
    EXPECT_NEAR(bl_recovered.r, bl_original.r, 1e-10);
    EXPECT_NEAR(bl_recovered.theta, bl_original.theta, 1e-10);
    EXPECT_NEAR(bl_recovered.phi, bl_original.phi, 1e-10);
}

// ============================================================================
// Numerical Precision Tests
// ============================================================================

TEST(PrecisionTest, NearHorizon) {
    // Metric should be well-defined very close to horizon
    KerrMetric<double> metric(1.0, 0.5);
    double r_horizon = metric.horizon_radius();
    double r_close = r_horizon * 1.001;
    
    EXPECT_NO_THROW({
        auto g = metric.compute(r_close, M_PI / 2.0);
        EXPECT_FALSE(std::isnan(g.g[0][0]));
        EXPECT_FALSE(std::isinf(g.g[1][1]));
    });
}

TEST(PrecisionTest, NearPoles) {
    // Metric should handle near-polar regions
    KerrMetric<double> metric(1.0, 0.5);
    
    EXPECT_NO_THROW({
        auto g1 = metric.compute(10.0, 0.001);
        auto g2 = metric.compute(10.0, M_PI - 0.001);
        EXPECT_FALSE(std::isnan(g1.g[3][3]));
        EXPECT_FALSE(std::isnan(g2.g[3][3]));
    });
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}