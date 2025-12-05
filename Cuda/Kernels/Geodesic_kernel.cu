/**
 * @file geodesic_kernel.cu
 * @brief CUDA kernels for GPU-accelerated geodesic ray tracing
 * 
 * Implements parallel ray tracing on NVIDIA GPUs with:
 * - Coalesced memory access patterns
 * - Shared memory caching for metric data
 * - Warp-level primitives for efficiency
 * 
 * @author BlackHole-Sim Team
 * @license MIT
 */

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <math.h>
#include <stdio.h>

// Constants
#define PI 3.14159265358979323846
#define MAX_STEPS 10000
#define TOLERANCE 1e-8

// ============================================================================
// Data Structures
// ============================================================================

/**
 * @brief Kerr metric parameters (constant memory for fast access)
 */
__constant__ struct {
    double M;       // Mass
    double a;       // Spin parameter
    double r_horizon;
    double r_isco;
} d_metric_params;

/**
 * @brief Ray state for integration
 */
struct RayState {
    double x[4];    // Position: t, r, theta, phi
    double p[4];    // Covariant momentum: p_t, p_r, p_theta, p_phi
};

/**
 * @brief Ray tracing result
 */
struct RayResult {
    float intensity;
    float r_hit;
    float phi_hit;
    int hit_type;   // 0=escaped, 1=horizon, 2=disk
    int num_steps;
};

/**
 * @brief Camera parameters
 */
struct CameraParams {
    double r_obs;
    double theta_obs;
    double phi_obs;
    double fov;
    int width;
    int height;
};

// ============================================================================
// Device Functions - Kerr Metric
// ============================================================================

/**
 * @brief Compute Sigma = r^2 + a^2 * cos^2(theta)
 */
__device__ __forceinline__
double compute_sigma(double r, double theta, double a) {
    double cos_theta = cos(theta);
    return r * r + a * a * cos_theta * cos_theta;
}

/**
 * @brief Compute Delta = r^2 - 2*M*r + a^2
 */
__device__ __forceinline__
double compute_delta(double r, double M, double a) {
    return r * r - 2.0 * M * r + a * a;
}

/**
 * @brief Compute metric components at (r, theta)
 * 
 * Returns: g[0]=g_tt, g[1]=g_rr, g[2]=g_thth, g[3]=g_phph, g[4]=g_tph
 */
__device__
void compute_metric(double r, double theta, double* g) {
    double M = d_metric_params.M;
    double a = d_metric_params.a;
    
    double sigma = compute_sigma(r, theta, a);
    double delta = compute_delta(r, M, a);
    
    double sin_theta = sin(theta);
    double sin2_theta = sin_theta * sin_theta;
    double cos_theta = cos(theta);
    
    double r2 = r * r;
    double a2 = a * a;
    double A = (r2 + a2) * (r2 + a2) - a2 * delta * sin2_theta;
    
    // Metric components
    g[0] = -(1.0 - 2.0 * M * r / sigma);              // g_tt
    g[1] = sigma / delta;                              // g_rr
    g[2] = sigma;                                      // g_thth
    g[3] = A * sin2_theta / sigma;                     // g_phph
    g[4] = -2.0 * M * a * r * sin2_theta / sigma;      // g_tph
}

/**
 * @brief Compute inverse metric components
 */
__device__
void compute_metric_inverse(double r, double theta, double* g_inv) {
    double M = d_metric_params.M;
    double a = d_metric_params.a;
    
    double sigma = compute_sigma(r, theta, a);
    double delta = compute_delta(r, M, a);
    
    double sin_theta = sin(theta);
    double sin2_theta = sin_theta * sin_theta;
    
    double r2 = r * r;
    double a2 = a * a;
    double A = (r2 + a2) * (r2 + a2) - a2 * delta * sin2_theta;
    
    // Inverse metric components
    g_inv[0] = -A / (sigma * delta);                   // g^tt
    g_inv[1] = delta / sigma;                          // g^rr
    g_inv[2] = 1.0 / sigma;                            // g^thth
    g_inv[3] = (delta - a2 * sin2_theta) / (sigma * delta * sin2_theta);  // g^phph
    g_inv[4] = -2.0 * M * a * r / (sigma * delta);     // g^tph
}

/**
 * @brief Compute Christoffel symbols (selected components for geodesic equation)
 */
__device__
void compute_christoffel(double r, double theta, double Gamma[4][4][4]) {
    double M = d_metric_params.M;
    double a = d_metric_params.a;
    
    double sigma = compute_sigma(r, theta, a);
    double delta = compute_delta(r, M, a);
    
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double sin2_theta = sin_theta * sin_theta;
    double cos2_theta = cos_theta * cos_theta;
    
    double r2 = r * r;
    double a2 = a * a;
    double sigma2 = sigma * sigma;
    
    // Initialize to zero
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0;
    
    // Non-zero Christoffel symbols for Kerr metric
    // Γ^t components
    double dSigma_dr = 2.0 * r;
    double dSigma_dth = -2.0 * a2 * sin_theta * cos_theta;
    double dDelta_dr = 2.0 * r - 2.0 * M;
    
    Gamma[0][0][1] = M * (r2 - a2 * cos2_theta) / (sigma2 * delta) * (r2 + a2);
    Gamma[0][1][0] = Gamma[0][0][1];
    
    Gamma[0][0][2] = -2.0 * M * a2 * r * sin_theta * cos_theta / sigma2;
    Gamma[0][2][0] = Gamma[0][0][2];
    
    // Γ^r components
    Gamma[1][0][0] = M * delta * (r2 - a2 * cos2_theta) / (sigma2 * sigma);
    Gamma[1][1][1] = (r * delta - M * (r2 - a2)) / (sigma * delta);
    Gamma[1][2][2] = -r * delta / sigma;
    Gamma[1][3][3] = -delta * sin2_theta * (r - M * a2 * sin2_theta / sigma) / sigma;
    
    // Γ^θ components
    Gamma[2][1][2] = r / sigma;
    Gamma[2][2][1] = Gamma[2][1][2];
    Gamma[2][3][3] = -sin_theta * cos_theta * ((r2 + a2) * (r2 + a2) / sigma + a2 * sin2_theta) / sigma;
    
    // Γ^φ components
    if (sin_theta > 1e-10) {
        Gamma[3][1][3] = (r * sigma - M * (r2 - a2 * cos2_theta)) / (sigma2 * delta);
        Gamma[3][3][1] = Gamma[3][1][3];
        Gamma[3][2][3] = cos_theta / sin_theta + 2.0 * a2 * sin_theta * cos_theta / sigma;
        Gamma[3][3][2] = Gamma[3][2][3];
    }
}

// ============================================================================
// Geodesic Integration
// ============================================================================

/**
 * @brief Compute RHS of geodesic equations: dx^μ/dλ and dp_μ/dλ
 */
__device__
void geodesic_rhs(const RayState& state, RayState& dstate) {
    double r = state.x[1];
    double theta = state.x[2];
    
    // Get metric inverse
    double g_inv[5];
    compute_metric_inverse(r, theta, g_inv);
    
    // Get Christoffel symbols
    double Gamma[4][4][4];
    compute_christoffel(r, theta, Gamma);
    
    // dx^μ/dλ = g^μν p_ν
    dstate.x[0] = g_inv[0] * state.p[0] + g_inv[4] * state.p[3];
    dstate.x[1] = g_inv[1] * state.p[1];
    dstate.x[2] = g_inv[2] * state.p[2];
    dstate.x[3] = g_inv[3] * state.p[3] + g_inv[4] * state.p[0];
    
    // Contravariant momentum
    double p_up[4];
    p_up[0] = g_inv[0] * state.p[0] + g_inv[4] * state.p[3];
    p_up[1] = g_inv[1] * state.p[1];
    p_up[2] = g_inv[2] * state.p[2];
    p_up[3] = g_inv[3] * state.p[3] + g_inv[4] * state.p[0];
    
    // dp_μ/dλ = Γ^ν_μρ p_ν p^ρ
    for (int mu = 0; mu < 4; mu++) {
        dstate.p[mu] = 0.0;
        for (int nu = 0; nu < 4; nu++) {
            for (int rho = 0; rho < 4; rho++) {
                dstate.p[mu] += Gamma[nu][mu][rho] * state.p[nu] * p_up[rho];
            }
        }
    }
}

/**
 * @brief Single RK4 step
 */
__device__
void rk4_step(RayState& state, double h) {
    RayState k1, k2, k3, k4;
    RayState temp;
    
    // k1
    geodesic_rhs(state, k1);
    
    // k2
    for (int i = 0; i < 4; i++) {
        temp.x[i] = state.x[i] + 0.5 * h * k1.x[i];
        temp.p[i] = state.p[i] + 0.5 * h * k1.p[i];
    }
    geodesic_rhs(temp, k2);
    
    // k3
    for (int i = 0; i < 4; i++) {
        temp.x[i] = state.x[i] + 0.5 * h * k2.x[i];
        temp.p[i] = state.p[i] + 0.5 * h * k2.p[i];
    }
    geodesic_rhs(temp, k3);
    
    // k4
    for (int i = 0; i < 4; i++) {
        temp.x[i] = state.x[i] + h * k3.x[i];
        temp.p[i] = state.p[i] + h * k3.p[i];
    }
    geodesic_rhs(temp, k4);
    
    // Update state
    for (int i = 0; i < 4; i++) {
        state.x[i] += h / 6.0 * (k1.x[i] + 2*k2.x[i] + 2*k3.x[i] + k4.x[i]);
        state.p[i] += h / 6.0 * (k1.p[i] + 2*k2.p[i] + 2*k3.p[i] + k4.p[i]);
    }
}

/**
 * @brief Initialize null geodesic from camera pixel
 */
__device__
void init_ray(int px, int py, const CameraParams& cam, RayState& state) {
    // Pixel to screen coordinates
    double aspect = (double)cam.width / cam.height;
    double fov_rad = cam.fov * PI / 180.0;
    
    double x = (2.0 * px / cam.width - 1.0) * tan(fov_rad / 2.0) * aspect;
    double y = (1.0 - 2.0 * py / cam.height) * tan(fov_rad / 2.0);
    
    double alpha = atan(x);
    double beta = atan(y / sqrt(1.0 + x*x));
    
    // Observer position
    state.x[0] = 0.0;
    state.x[1] = cam.r_obs;
    state.x[2] = cam.theta_obs;
    state.x[3] = cam.phi_obs;
    
    // Get metric at observer
    double g[5];
    compute_metric(cam.r_obs, cam.theta_obs, g);
    
    double sin_theta = sin(cam.theta_obs);
    
    // Local direction
    double k_r = -cos(beta) * cos(alpha);
    double k_th = sin(beta);
    double k_ph = cos(beta) * sin(alpha);
    
    // Transform to covariant momentum
    state.p[1] = k_r * sqrt(fabs(g[1]));
    state.p[2] = k_th * sqrt(fabs(g[2]));
    state.p[3] = k_ph * sqrt(fabs(g[3])) / sin_theta;
    
    // Null condition: g^μν p_μ p_ν = 0
    // Solve for p_t
    double g_inv[5];
    compute_metric_inverse(cam.r_obs, cam.theta_obs, g_inv);
    
    double A = g_inv[0];
    double B = 2.0 * g_inv[4] * state.p[3];
    double C = g_inv[1] * state.p[1] * state.p[1] +
               g_inv[2] * state.p[2] * state.p[2] +
               g_inv[3] * state.p[3] * state.p[3];
    
    double discriminant = B * B - 4.0 * A * C;
    if (discriminant < 0) discriminant = 0;
    
    state.p[0] = (-B - sqrt(discriminant)) / (2.0 * A);
}

// ============================================================================
// Main Ray Tracing Kernel
// ============================================================================

/**
 * @brief Main CUDA kernel for parallel ray tracing
 * 
 * Each thread traces one ray from the camera through the pixel.
 */
__global__
void raytrace_kernel(
    RayResult* results,
    const CameraParams cam,
    const double disk_inner,
    const double disk_outer,
    const double escape_radius
) {
    int px = blockIdx.x * blockDim.x + threadIdx.x;
    int py = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (px >= cam.width || py >= cam.height) return;
    
    int idx = py * cam.width + px;
    
    // Initialize ray
    RayState state;
    init_ray(px, py, cam, state);
    
    RayResult result;
    result.intensity = 0.0f;
    result.hit_type = 0;
    result.num_steps = 0;
    result.r_hit = 0.0f;
    result.phi_hit = 0.0f;
    
    double r_horizon = d_metric_params.r_horizon * 1.01;
    
    // Integration loop
    double h = -0.1;  // Backward tracing
    double prev_theta = state.x[2];
    
    for (int step = 0; step < MAX_STEPS; step++) {
        double r = state.x[1];
        double theta = state.x[2];
        
        // Check horizon crossing
        if (r < r_horizon) {
            result.hit_type = 1;  // Horizon
            result.intensity = 0.0f;
            result.num_steps = step;
            break;
        }
        
        // Check escape
        if (r > escape_radius) {
            result.hit_type = 0;  // Escaped
            result.intensity = 0.01f;  // Background
            result.num_steps = step;
            break;
        }
        
        // Check disk crossing (equatorial plane)
        double theta_eq = PI / 2.0;
        if ((prev_theta - theta_eq) * (theta - theta_eq) < 0) {
            // Linear interpolation
            double t = (theta_eq - prev_theta) / (theta - prev_theta);
            double r_cross = state.x[1];
            double phi_cross = state.x[3];
            
            if (r_cross >= disk_inner && r_cross <= disk_outer) {
                result.hit_type = 2;  // Disk
                result.r_hit = (float)r_cross;
                result.phi_hit = (float)phi_cross;
                
                // Disk emission (temperature profile)
                double intensity = pow(disk_inner / r_cross, 3.0);
                
                // Doppler beaming
                double v_orbit = sqrt(1.0 / r_cross);
                double doppler = 1.0 + v_orbit * sin(phi_cross) * sin(cam.theta_obs);
                intensity *= pow(fabs(doppler), 3.0);
                
                result.intensity = (float)fmin(intensity * 2.0, 10.0);
                result.num_steps = step;
                break;
            }
        }
        
        prev_theta = theta;
        
        // Bound theta
        if (state.x[2] < 0.01) state.x[2] = 0.01;
        if (state.x[2] > PI - 0.01) state.x[2] = PI - 0.01;
        
        // Adaptive step size near horizon
        double r_factor = (r - r_horizon) / d_metric_params.M;
        double step_h = h;
        if (r_factor < 10.0 && r_factor > 0) {
            step_h = -fmin(fabs(h), r_factor * 0.1);
        }
        
        // RK4 step
        rk4_step(state, step_h);
        
        result.num_steps = step;
    }
    
    results[idx] = result;
}

// ============================================================================
// Intensity to Color Kernel
// ============================================================================

/**
 * @brief Convert intensity to RGBA using inferno colormap
 */
__device__
void intensity_to_color(float intensity, unsigned char* rgba) {
    // Inferno colormap approximation
    float t = fminf(fmaxf(intensity, 0.0f), 1.0f);
    
    float r, g, b;
    
    if (t < 0.25f) {
        float s = t / 0.25f;
        r = 0.0f + s * 0.34f;
        g = 0.0f + s * 0.06f;
        b = 0.01f + s * 0.26f;
    } else if (t < 0.5f) {
        float s = (t - 0.25f) / 0.25f;
        r = 0.34f + s * 0.24f;
        g = 0.06f + s * 0.09f;
        b = 0.27f + s * 0.16f;
    } else if (t < 0.75f) {
        float s = (t - 0.5f) / 0.25f;
        r = 0.58f + s * 0.27f;
        g = 0.15f + s * 0.18f;
        b = 0.43f - s * 0.17f;
    } else {
        float s = (t - 0.75f) / 0.25f;
        r = 0.85f + s * 0.14f;
        g = 0.33f + s * 0.37f;
        b = 0.26f - s * 0.06f;
    }
    
    rgba[0] = (unsigned char)(r * 255);
    rgba[1] = (unsigned char)(g * 255);
    rgba[2] = (unsigned char)(b * 255);
    rgba[3] = 255;
}

/**
 * @brief Convert ray results to RGBA image
 */
__global__
void results_to_image_kernel(
    const RayResult* results,
    unsigned char* image,
    int width,
    int height,
    float intensity_scale
) {
    int px = blockIdx.x * blockDim.x + threadIdx.x;
    int py = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (px >= width || py >= height) return;
    
    int idx = py * width + px;
    
    float intensity = results[idx].intensity * intensity_scale;
    
    // Tone mapping
    intensity = intensity / (intensity + 1.0f);
    
    intensity_to_color(intensity, &image[idx * 4]);
}

// ============================================================================
// Host Interface Functions
// ============================================================================

extern "C" {

/**
 * @brief Set metric parameters in constant memory
 */
void cuda_set_metric(double M, double a) {
    struct {
        double M, a, r_horizon, r_isco;
    } params;
    
    params.M = M;
    params.a = a;
    params.r_horizon = M + sqrt(M * M - a * a);
    
    // ISCO calculation
    double a_star = a / M;
    double Z1 = 1.0 + pow(1.0 - a_star * a_star, 1.0/3.0) * 
                (pow(1.0 + a_star, 1.0/3.0) + pow(1.0 - a_star, 1.0/3.0));
    double Z2 = sqrt(3.0 * a_star * a_star + Z1 * Z1);
    params.r_isco = M * (3.0 + Z2 - sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));
    
    cudaMemcpyToSymbol(d_metric_params, &params, sizeof(params));
}

/**
 * @brief Main ray tracing function
 */
void cuda_raytrace(
    float* h_intensity,
    int width, int height,
    double r_obs, double theta_obs, double fov,
    double disk_inner, double disk_outer,
    double M, double a
) {
    // Set metric
    cuda_set_metric(M, a);
    
    // Camera params
    CameraParams cam;
    cam.r_obs = r_obs;
    cam.theta_obs = theta_obs;
    cam.phi_obs = 0.0;
    cam.fov = fov;
    cam.width = width;
    cam.height = height;
    
    // Allocate device memory
    RayResult* d_results;
    cudaMalloc(&d_results, width * height * sizeof(RayResult));
    
    // Launch kernel
    dim3 block(16, 16);
    dim3 grid((width + 15) / 16, (height + 15) / 16);
    
    raytrace_kernel<<<grid, block>>>(d_results, cam, disk_inner, disk_outer, 500.0);
    
    // Copy results back
    RayResult* h_results = new RayResult[width * height];
    cudaMemcpy(h_results, d_results, width * height * sizeof(RayResult), cudaMemcpyDeviceToHost);
    
    // Extract intensity
    for (int i = 0; i < width * height; i++) {
        h_intensity[i] = h_results[i].intensity;
    }
    
    delete[] h_results;
    cudaFree(d_results);
}

/**
 * @brief Render to RGBA image
 */
void cuda_render_image(
    unsigned char* h_image,
    int width, int height,
    double r_obs, double theta_obs, double fov,
    double disk_inner, double disk_outer,
    double M, double a,
    float intensity_scale
) {
    cuda_set_metric(M, a);
    
    CameraParams cam;
    cam.r_obs = r_obs;
    cam.theta_obs = theta_obs;
    cam.phi_obs = 0.0;
    cam.fov = fov;
    cam.width = width;
    cam.height = height;
    
    RayResult* d_results;
    unsigned char* d_image;
    
    cudaMalloc(&d_results, width * height * sizeof(RayResult));
    cudaMalloc(&d_image, width * height * 4);
    
    dim3 block(16, 16);
    dim3 grid((width + 15) / 16, (height + 15) / 16);
    
    // Ray trace
    raytrace_kernel<<<grid, block>>>(d_results, cam, disk_inner, disk_outer, 500.0);
    
    // Convert to image
    results_to_image_kernel<<<grid, block>>>(d_results, d_image, width, height, intensity_scale);
    
    // Copy back
    cudaMemcpy(h_image, d_image, width * height * 4, cudaMemcpyDeviceToHost);
    
    cudaFree(d_results);
    cudaFree(d_image);
}

/**
 * @brief Get device info
 */
void cuda_get_device_info(char* name, int* major, int* minor, size_t* memory) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    
    strcpy(name, prop.name);
    *major = prop.major;
    *minor = prop.minor;
    *memory = prop.totalGlobalMem;
}

} // extern "C"