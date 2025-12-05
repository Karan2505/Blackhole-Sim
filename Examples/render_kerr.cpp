/**
 * @file render_kerr.cpp
 * @brief Example: Render a Kerr black hole image
 * 
 * This example demonstrates:
 * - Creating a Kerr metric with custom spin
 * - Setting up a camera for ray tracing
 * - Rendering an image with accretion disk
 * - Saving output as PPM image
 * 
 * Usage:
 *   ./render_kerr [spin] [resolution] [output.ppm]
 *   ./render_kerr 0.99 512 gargantua.ppm
 * 
 * @author BlackHole-Sim Team
 * @license MIT
 */

#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cmath>

#include "blackhole/metric.hpp"
#include "blackhole/geodesic.hpp"

using namespace blackhole;

/**
 * @brief Inferno colormap (scientific visualization)
 */
void inferno_colormap(double t, unsigned char& r, unsigned char& g, unsigned char& b) {
    t = std::max(0.0, std::min(1.0, t));
    
    // Inferno colormap coefficients
    double r_val, g_val, b_val;
    
    if (t < 0.25) {
        double s = t / 0.25;
        r_val = 0.0 + s * 0.34;
        g_val = 0.0 + s * 0.06;
        b_val = 0.01 + s * 0.26;
    } else if (t < 0.5) {
        double s = (t - 0.25) / 0.25;
        r_val = 0.34 + s * 0.24;
        g_val = 0.06 + s * 0.09;
        b_val = 0.27 + s * 0.16;
    } else if (t < 0.75) {
        double s = (t - 0.5) / 0.25;
        r_val = 0.58 + s * 0.27;
        g_val = 0.15 + s * 0.18;
        b_val = 0.43 - s * 0.17;
    } else {
        double s = (t - 0.75) / 0.25;
        r_val = 0.85 + s * 0.14;
        g_val = 0.33 + s * 0.37;
        b_val = 0.26 - s * 0.06;
    }
    
    r = static_cast<unsigned char>(r_val * 255);
    g = static_cast<unsigned char>(g_val * 255);
    b = static_cast<unsigned char>(b_val * 255);
}

/**
 * @brief Save image as PPM file
 */
void save_ppm(const std::string& filename, 
              const std::vector<std::vector<double>>& image,
              int width, int height) {
    std::ofstream file(filename, std::ios::binary);
    file << "P6\n" << width << " " << height << "\n255\n";
    
    // Find max intensity for normalization
    double max_intensity = 0;
    for (const auto& row : image) {
        for (double val : row) {
            max_intensity = std::max(max_intensity, val);
        }
    }
    if (max_intensity < 1e-10) max_intensity = 1.0;
    
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            double intensity = image[j][i] / max_intensity;
            
            // Tone mapping
            intensity = intensity / (intensity + 1.0);
            
            // Gamma correction
            intensity = std::pow(intensity, 1.0 / 2.2);
            
            unsigned char r, g, b;
            inferno_colormap(intensity, r, g, b);
            
            file.put(r);
            file.put(g);
            file.put(b);
        }
    }
    
    file.close();
}

int main(int argc, char* argv[]) {
    // Parse arguments
    double spin = 0.99;
    int resolution = 512;
    std::string output = "kerr_shadow.ppm";
    
    if (argc >= 2) spin = std::atof(argv[1]);
    if (argc >= 3) resolution = std::atoi(argv[2]);
    if (argc >= 4) output = argv[3];
    
    // Validate spin
    if (std::abs(spin) >= 1.0) {
        std::cerr << "Error: Spin must be |a/M| < 1. Got: " << spin << std::endl;
        return 1;
    }
    
    std::cout << "╔══════════════════════════════════════════════╗" << std::endl;
    std::cout << "║     BlackHole-Sim: Kerr Black Hole Render    ║" << std::endl;
    std::cout << "╚══════════════════════════════════════════════╝" << std::endl;
    std::cout << std::endl;
    
    // Create Kerr metric
    KerrMetric<double> metric(1.0, spin);
    
    // Print black hole properties
    std::cout << "Black Hole Properties:" << std::endl;
    std::cout << "  Mass (M):           1.0 (geometric units)" << std::endl;
    std::cout << "  Spin (a/M):         " << spin << std::endl;
    std::cout << "  Event Horizon:      " << metric.horizon_radius() << " M" << std::endl;
    std::cout << "  ISCO (prograde):    " << metric.isco_radius() << " M" << std::endl;
    std::cout << "  ISCO (retrograde):  " << metric.isco_radius_retrograde() << " M" << std::endl;
    std::cout << "  Photon Sphere:      " << metric.photon_sphere_radius() << " M" << std::endl;
    std::cout << "  Ω_H (horizon):      " << metric.horizon_angular_velocity() << " c/M" << std::endl;
    std::cout << std::endl;
    
    // Camera setup
    double r_observer = 100.0;     // Observer distance
    double theta_observer = 85.0 * M_PI / 180.0;  // Near edge-on (like Interstellar)
    double fov = 15.0;             // Field of view in degrees
    
    Camera<double> camera(r_observer, theta_observer, fov, resolution, resolution);
    
    std::cout << "Camera Setup:" << std::endl;
    std::cout << "  Resolution:    " << resolution << " x " << resolution << std::endl;
    std::cout << "  Observer r:    " << r_observer << " M" << std::endl;
    std::cout << "  Inclination:   " << 90.0 - theta_observer * 180.0 / M_PI << "°" << std::endl;
    std::cout << "  Field of View: " << fov << "°" << std::endl;
    std::cout << std::endl;
    
    // Accretion disk parameters
    double disk_inner = metric.isco_radius();
    double disk_outer = 20.0;
    
    std::cout << "Accretion Disk:" << std::endl;
    std::cout << "  Inner Edge: " << disk_inner << " M (at ISCO)" << std::endl;
    std::cout << "  Outer Edge: " << disk_outer << " M" << std::endl;
    std::cout << std::endl;
    
    // Integrator configuration
    IntegratorConfig<double> config;
    config.max_affine_param = 500.0;
    config.tolerance = 1e-8;
    config.store_trajectory = false;  // Save memory
    
    // Ray tracer
    RayTracer<double> tracer(metric, camera, config);
    
    std::cout << "Rendering..." << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Render image
    auto image = tracer.render(disk_inner, disk_outer);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    // Statistics
    int total_rays = resolution * resolution;
    double rays_per_second = total_rays / (duration.count() / 1000.0);
    
    std::cout << std::endl;
    std::cout << "Render Statistics:" << std::endl;
    std::cout << "  Total Rays:     " << total_rays << std::endl;
    std::cout << "  Render Time:    " << duration.count() / 1000.0 << " seconds" << std::endl;
    std::cout << "  Rays/Second:    " << rays_per_second << std::endl;
    std::cout << std::endl;
    
    // Save image
    save_ppm(output, image, resolution, resolution);
    
    std::cout << "Output saved to: " << output << std::endl;
    std::cout << std::endl;
    std::cout << "Done! 🌌" << std::endl;
    
    return 0;
}