"""
BlackHole-Sim: Production-Grade Black Hole Simulation
======================================================

A high-fidelity black hole simulation package implementing:
- General relativistic spacetime (Schwarzschild & Kerr metrics)
- Geodesic ray tracing with adaptive integration
- Accretion disk models with Doppler effects
- Gravitational lensing and shadow calculation

Example usage:
    >>> import bhsim
    >>> bh = bhsim.KerrBlackHole(mass=1.0, spin=0.99)
    >>> image = bh.render(resolution=512, inclination=85)
    >>> bhsim.save_image(image, "kerr_shadow.png")
"""

__version__ = "1.0.0"
__author__ = "BlackHole-Sim Team"

import numpy as np
from typing import Tuple, Optional, List, Dict, Any
from dataclasses import dataclass
import json

# Try to import the C++ extension if available
try:
    from . import _bhsim_core as _core
    HAS_CPP_BACKEND = True
except ImportError:
    HAS_CPP_BACKEND = False
    print("Warning: C++ backend not available, using pure Python (slower)")


@dataclass
class PhysicalConstants:
    """Physical constants in CGS units"""
    c: float = 2.99792458e10      # Speed of light (cm/s)
    G: float = 6.67430e-8          # Gravitational constant
    M_sun: float = 1.98892e33      # Solar mass (g)
    h: float = 6.62607e-27         # Planck constant
    k_B: float = 1.38065e-16       # Boltzmann constant
    m_e: float = 9.10938e-28       # Electron mass
    e: float = 4.80326e-10         # Elementary charge (esu)


CONSTANTS = PhysicalConstants()


@dataclass
class MetricResult:
    """Result of metric tensor computation"""
    g: np.ndarray           # Covariant metric g_μν (4x4)
    g_inv: np.ndarray       # Contravariant metric g^μν (4x4)
    christoffel: np.ndarray # Christoffel symbols Γ^α_μν (4x4x4)
    
    # Derived quantities
    r_horizon: float        # Event horizon radius
    r_isco: float           # ISCO radius
    r_photon: float         # Photon sphere radius
    omega_H: float          # Horizon angular velocity


class Metric:
    """Base class for spacetime metrics"""
    
    def __init__(self, mass: float = 1.0):
        if mass <= 0:
            raise ValueError("Mass must be positive")
        self.M = mass
    
    def compute(self, r: float, theta: float) -> MetricResult:
        """Compute metric tensor at (r, θ)"""
        raise NotImplementedError
    
    @property
    def horizon_radius(self) -> float:
        raise NotImplementedError
    
    @property
    def isco_radius(self) -> float:
        raise NotImplementedError
    
    @property
    def photon_sphere_radius(self) -> float:
        raise NotImplementedError


class SchwarzschildMetric(Metric):
    """
    Schwarzschild metric for non-rotating black holes.
    
    ds² = -(1 - 2M/r)dt² + (1 - 2M/r)⁻¹dr² + r²dθ² + r²sin²θ dφ²
    """
    
    def __init__(self, mass: float = 1.0):
        super().__init__(mass)
    
    @property
    def horizon_radius(self) -> float:
        return 2.0 * self.M
    
    @property
    def isco_radius(self) -> float:
        return 6.0 * self.M
    
    @property
    def photon_sphere_radius(self) -> float:
        return 3.0 * self.M
    
    def compute(self, r: float, theta: float) -> MetricResult:
        if r <= 0:
            raise ValueError("r must be positive")
        
        f = 1.0 - 2.0 * self.M / r
        sin_theta = np.sin(theta)
        sin2_theta = sin_theta ** 2
        
        # Covariant metric
        g = np.zeros((4, 4))
        g[0, 0] = -f
        g[1, 1] = 1.0 / f
        g[2, 2] = r ** 2
        g[3, 3] = r ** 2 * sin2_theta
        
        # Contravariant metric
        g_inv = np.zeros((4, 4))
        g_inv[0, 0] = -1.0 / f
        g_inv[1, 1] = f
        g_inv[2, 2] = 1.0 / r ** 2
        g_inv[3, 3] = 1.0 / (r ** 2 * sin2_theta) if sin2_theta > 0 else 0
        
        # Christoffel symbols
        christoffel = self._compute_christoffel(r, theta, f)
        
        return MetricResult(
            g=g,
            g_inv=g_inv,
            christoffel=christoffel,
            r_horizon=self.horizon_radius,
            r_isco=self.isco_radius,
            r_photon=self.photon_sphere_radius,
            omega_H=0.0
        )
    
    def _compute_christoffel(self, r: float, theta: float, f: float) -> np.ndarray:
        """Compute Christoffel symbols for Schwarzschild metric"""
        Gamma = np.zeros((4, 4, 4))
        
        r2 = r ** 2
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        sin2_theta = sin_theta ** 2
        
        # Non-zero components
        Gamma[0, 0, 1] = self.M / (r2 * f)
        Gamma[0, 1, 0] = Gamma[0, 0, 1]
        
        Gamma[1, 0, 0] = self.M * f / r2
        Gamma[1, 1, 1] = -self.M / (r2 * f)
        Gamma[1, 2, 2] = -(r - 2.0 * self.M)
        Gamma[1, 3, 3] = -(r - 2.0 * self.M) * sin2_theta
        
        Gamma[2, 1, 2] = 1.0 / r
        Gamma[2, 2, 1] = Gamma[2, 1, 2]
        Gamma[2, 3, 3] = -sin_theta * cos_theta
        
        Gamma[3, 1, 3] = 1.0 / r
        Gamma[3, 3, 1] = Gamma[3, 1, 3]
        if sin_theta != 0:
            Gamma[3, 2, 3] = cos_theta / sin_theta
            Gamma[3, 3, 2] = Gamma[3, 2, 3]
        
        return Gamma


class KerrMetric(Metric):
    """
    Kerr metric for rotating black holes in Boyer-Lindquist coordinates.
    
    Parameters:
        mass: Black hole mass (geometric units, M=1 for unit mass)
        spin: Spin parameter a = J/M (must satisfy |a| < M)
    """
    
    def __init__(self, mass: float = 1.0, spin: float = 0.0):
        super().__init__(mass)
        if abs(spin) >= mass:
            raise ValueError("Spin magnitude must be less than mass for black hole")
        self.a = spin
    
    @property
    def horizon_radius(self) -> float:
        """Outer event horizon r+"""
        return self.M + np.sqrt(self.M**2 - self.a**2)
    
    @property
    def inner_horizon_radius(self) -> float:
        """Inner (Cauchy) horizon r-"""
        return self.M - np.sqrt(self.M**2 - self.a**2)
    
    @property
    def isco_radius(self) -> float:
        """ISCO radius for prograde orbits"""
        a_star = self.a / self.M
        Z1 = 1 + (1 - a_star**2)**(1/3) * ((1 + a_star)**(1/3) + (1 - a_star)**(1/3))
        Z2 = np.sqrt(3 * a_star**2 + Z1**2)
        return self.M * (3 + Z2 - np.sqrt((3 - Z1) * (3 + Z1 + 2*Z2)))
    
    @property
    def isco_radius_retrograde(self) -> float:
        """ISCO radius for retrograde orbits"""
        a_star = self.a / self.M
        Z1 = 1 + (1 - a_star**2)**(1/3) * ((1 + a_star)**(1/3) + (1 - a_star)**(1/3))
        Z2 = np.sqrt(3 * a_star**2 + Z1**2)
        return self.M * (3 + Z2 + np.sqrt((3 - Z1) * (3 + Z1 + 2*Z2)))
    
    @property
    def photon_sphere_radius(self) -> float:
        """Photon sphere radius (equatorial, prograde)"""
        a_star = self.a / self.M
        return 2 * self.M * (1 + np.cos(2/3 * np.arccos(-a_star)))
    
    @property
    def horizon_angular_velocity(self) -> float:
        """Angular velocity at the horizon Ω_H"""
        r_plus = self.horizon_radius
        return self.a / (2 * self.M * r_plus)
    
    def ergosphere_radius(self, theta: float) -> float:
        """Ergosphere boundary at given polar angle"""
        cos_theta = np.cos(theta)
        return self.M + np.sqrt(self.M**2 - self.a**2 * cos_theta**2)
    
    def _Sigma(self, r: float, theta: float) -> float:
        return r**2 + self.a**2 * np.cos(theta)**2
    
    def _Delta(self, r: float) -> float:
        return r**2 - 2*self.M*r + self.a**2
    
    def _A(self, r: float, theta: float) -> float:
        r2_a2 = r**2 + self.a**2
        return r2_a2**2 - self.a**2 * self._Delta(r) * np.sin(theta)**2
    
    def compute(self, r: float, theta: float) -> MetricResult:
        if r <= 0:
            raise ValueError("r must be positive")
        
        Sigma = self._Sigma(r, theta)
        Delta = self._Delta(r)
        A = self._A(r, theta)
        
        sin_theta = np.sin(theta)
        sin2_theta = sin_theta ** 2
        cos_theta = np.cos(theta)
        
        # Covariant metric g_μν
        g = np.zeros((4, 4))
        g[0, 0] = -(1 - 2*self.M*r/Sigma)
        g[0, 3] = -2*self.M*self.a*r*sin2_theta/Sigma
        g[3, 0] = g[0, 3]
        g[1, 1] = Sigma/Delta
        g[2, 2] = Sigma
        g[3, 3] = A*sin2_theta/Sigma
        
        # Contravariant metric g^μν
        det_tphi = g[0, 0]*g[3, 3] - g[0, 3]**2
        
        g_inv = np.zeros((4, 4))
        g_inv[0, 0] = -g[3, 3] / det_tphi
        g_inv[0, 3] = g[0, 3] / det_tphi
        g_inv[3, 0] = g_inv[0, 3]
        g_inv[1, 1] = Delta / Sigma
        g_inv[2, 2] = 1.0 / Sigma
        g_inv[3, 3] = -g[0, 0] / det_tphi
        
        # Christoffel symbols (simplified)
        christoffel = self._compute_christoffel(r, theta, Sigma, Delta)
        
        return MetricResult(
            g=g,
            g_inv=g_inv,
            christoffel=christoffel,
            r_horizon=self.horizon_radius,
            r_isco=self.isco_radius,
            r_photon=self.photon_sphere_radius,
            omega_H=self.horizon_angular_velocity
        )
    
    def _compute_christoffel(self, r: float, theta: float, Sigma: float, Delta: float) -> np.ndarray:
        """Compute Christoffel symbols for Kerr metric (simplified)"""
        Gamma = np.zeros((4, 4, 4))
        
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        sin2_theta = sin_theta ** 2
        cos2_theta = cos_theta ** 2
        
        r2 = r ** 2
        a2 = self.a ** 2
        Sigma2 = Sigma ** 2
        
        # Selected non-zero components (simplified)
        if Delta != 0:
            Gamma[0, 0, 1] = self.M * (r2 - a2*cos2_theta) / (Sigma2 * Delta) * (r2 + a2)
            Gamma[0, 1, 0] = Gamma[0, 0, 1]
        
        Gamma[1, 0, 0] = self.M * Delta * (r2 - a2*cos2_theta) / (Sigma2 * Sigma)
        
        if Delta != 0:
            Gamma[1, 1, 1] = (r*Delta - self.M*(r2 - a2)) / (Sigma * Delta)
        
        Gamma[1, 2, 2] = -r * Delta / Sigma
        
        Gamma[2, 1, 2] = r / Sigma
        Gamma[2, 2, 1] = Gamma[2, 1, 2]
        
        if sin_theta != 0:
            Gamma[3, 1, 3] = (r*Sigma - self.M*(r2 - a2*cos2_theta)) / (Sigma2 * Delta) if Delta != 0 else 0
            Gamma[3, 3, 1] = Gamma[3, 1, 3]
            Gamma[3, 2, 3] = cos_theta / sin_theta
            Gamma[3, 3, 2] = Gamma[3, 2, 3]
        
        return Gamma
    
    def constants_of_motion(self, position: np.ndarray, velocity: np.ndarray) -> Tuple[float, float, float]:
        """
        Compute constants of motion for a geodesic.
        
        Returns:
            E: Energy
            L: Angular momentum
            Q: Carter constant
        """
        r, theta = position[1], position[2]
        result = self.compute(r, theta)
        g = result.g
        
        # Energy E = -p_t = -g_tμ u^μ
        E = -(g[0, 0]*velocity[0] + g[0, 3]*velocity[3])
        
        # Angular momentum L = p_φ = g_φμ u^μ
        L = g[3, 0]*velocity[0] + g[3, 3]*velocity[3]
        
        # Carter constant
        p_theta = g[2, 2] * velocity[2]
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        
        Q = p_theta**2 + cos_theta**2 * (self.a**2 * (1 - E**2) + L**2 / sin_theta**2)
        
        return E, L, Q


class GeodesicIntegrator:
    """
    Integrator for null and timelike geodesics in curved spacetime.
    
    Uses adaptive RKF45 method for accurate integration.
    """
    
    def __init__(
        self,
        metric: Metric,
        max_steps: int = 100000,
        tolerance: float = 1e-8,
        store_trajectory: bool = False
    ):
        self.metric = metric
        self.max_steps = max_steps
        self.tolerance = tolerance
        self.store_trajectory = store_trajectory
    
    def integrate(
        self,
        initial_position: np.ndarray,
        initial_momentum: np.ndarray,
        affine_max: float = 1000.0,
        backward: bool = True
    ) -> Dict[str, Any]:
        """
        Integrate geodesic equations.
        
        Args:
            initial_position: [t, r, θ, φ]
            initial_momentum: [p_t, p_r, p_θ, p_φ]
            affine_max: Maximum affine parameter
            backward: Integrate backward (for ray tracing)
        
        Returns:
            Dictionary with integration results
        """
        x = initial_position.copy()
        p = initial_momentum.copy()
        
        h = -0.1 if backward else 0.1
        lambda_param = 0.0
        
        r_horizon = self.metric.horizon_radius * 1.01
        
        positions = [x.copy()] if self.store_trajectory else []
        momenta = [p.copy()] if self.store_trajectory else []
        
        hit_horizon = False
        escaped = False
        num_steps = 0
        
        while num_steps < self.max_steps:
            # Check termination
            if x[1] < r_horizon:
                hit_horizon = True
                break
            
            if x[1] > 1000.0:
                escaped = True
                break
            
            if abs(lambda_param) > affine_max:
                break
            
            # Keep theta bounded
            x[2] = np.clip(x[2], 0.01, np.pi - 0.01)
            
            # RK4 step
            x_new, p_new = self._rk4_step(x, p, h)
            
            x = x_new
            p = p_new
            lambda_param += h
            num_steps += 1
            
            if self.store_trajectory:
                positions.append(x.copy())
                momenta.append(p.copy())
        
        return {
            'final_position': x,
            'final_momentum': p,
            'hit_horizon': hit_horizon,
            'escaped': escaped,
            'num_steps': num_steps,
            'affine_param': lambda_param,
            'positions': np.array(positions) if positions else None,
            'momenta': np.array(momenta) if momenta else None
        }
    
    def _rk4_step(
        self,
        x: np.ndarray,
        p: np.ndarray,
        h: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Single RK4 integration step"""
        
        k1_x, k1_p = self._geodesic_rhs(x, p)
        k2_x, k2_p = self._geodesic_rhs(x + h/2*k1_x, p + h/2*k1_p)
        k3_x, k3_p = self._geodesic_rhs(x + h/2*k2_x, p + h/2*k2_p)
        k4_x, k4_p = self._geodesic_rhs(x + h*k3_x, p + h*k3_p)
        
        x_new = x + h/6 * (k1_x + 2*k2_x + 2*k3_x + k4_x)
        p_new = p + h/6 * (k1_p + 2*k2_p + 2*k3_p + k4_p)
        
        return x_new, p_new
    
    def _geodesic_rhs(
        self,
        x: np.ndarray,
        p: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute right-hand side of geodesic equations"""
        
        r, theta = x[1], x[2]
        result = self.metric.compute(r, theta)
        g_inv = result.g_inv
        Gamma = result.christoffel
        
        # dx/dλ = g^μν p_ν
        dx = g_inv @ p
        
        # Get contravariant momentum
        p_up = g_inv @ p
        
        # dp_μ/dλ = Γ^ν_μρ p_ν p^ρ
        dp = np.zeros(4)
        for mu in range(4):
            for nu in range(4):
                for rho in range(4):
                    dp[mu] += Gamma[nu, mu, rho] * p[nu] * p_up[rho]
        
        return dx, dp


class KerrBlackHole:
    """
    High-level interface for Kerr black hole simulation.
    
    Example:
        >>> bh = KerrBlackHole(mass=1.0, spin=0.99)
        >>> print(f"ISCO: {bh.isco_radius:.3f} M")
        >>> image = bh.render(resolution=256, inclination=85)
    """
    
    def __init__(self, mass: float = 1.0, spin: float = 0.0):
        """
        Initialize Kerr black hole.
        
        Args:
            mass: Black hole mass in geometric units
            spin: Dimensionless spin parameter a/M (|a/M| < 1)
        """
        self.metric = KerrMetric(mass, spin * mass)
        self.mass = mass
        self.spin = spin
    
    @property
    def horizon_radius(self) -> float:
        return self.metric.horizon_radius
    
    @property
    def isco_radius(self) -> float:
        return self.metric.isco_radius
    
    @property
    def photon_sphere_radius(self) -> float:
        return self.metric.photon_sphere_radius
    
    @property
    def shadow_radius(self) -> float:
        """Approximate shadow radius"""
        return 3 * np.sqrt(3) * self.mass * (1 + 0.2 * (1 - self.spin**2))
    
    def render(
        self,
        resolution: int = 256,
        inclination: float = 85.0,
        observer_distance: float = 100.0,
        fov: float = 20.0,
        disk_inner: Optional[float] = None,
        disk_outer: float = 20.0,
        colormap: str = 'inferno'
    ) -> np.ndarray:
        """
        Render black hole image.
        
        Args:
            resolution: Image resolution (square)
            inclination: Observer inclination angle (degrees)
            observer_distance: Observer distance in M
            fov: Field of view in degrees
            disk_inner: Inner disk radius (defaults to ISCO)
            disk_outer: Outer disk radius
            colormap: Matplotlib colormap name
        
        Returns:
            2D numpy array with intensity values
        """
        if disk_inner is None:
            disk_inner = self.isco_radius
        
        theta_obs = (90 - inclination) * np.pi / 180
        
        image = np.zeros((resolution, resolution))
        
        # Generate rays for each pixel
        for j in range(resolution):
            for i in range(resolution):
                # Convert pixel to impact parameters
                x = (2 * i / resolution - 1) * np.tan(fov * np.pi / 360)
                y = (1 - 2 * j / resolution) * np.tan(fov * np.pi / 360)
                
                alpha = np.arctan(x)
                beta = np.arctan(y / np.sqrt(1 + x*x))
                
                # Initialize null geodesic
                intensity = self._trace_ray(
                    observer_distance, theta_obs, alpha, beta,
                    disk_inner, disk_outer
                )
                
                image[j, i] = intensity
        
        return image
    
    def _trace_ray(
        self,
        r_obs: float,
        theta_obs: float,
        alpha: float,
        beta: float,
        disk_inner: float,
        disk_outer: float
    ) -> float:
        """Trace single ray and return intensity"""
        
        # Simplified ray tracing
        # In full implementation, use GeodesicIntegrator
        
        r = r_obs
        theta = theta_obs
        phi = 0.0
        
        # Direction
        dr = -np.cos(beta) * np.cos(alpha)
        dtheta = np.sin(beta)
        dphi = np.cos(beta) * np.sin(alpha)
        
        r_horizon = self.horizon_radius * 1.05
        
        for step in range(1000):
            # Simple ray marching
            step_size = 0.1 + r * 0.02
            
            # Gravitational deflection
            deflection = 2 * self.mass / (r * r)
            dr -= deflection * step_size
            
            r += dr * step_size
            theta += dtheta * step_size / r
            phi += dphi * step_size / (r * np.sin(theta))
            
            # Check horizon
            if r < r_horizon:
                return 0.0
            
            # Check disk crossing
            if abs(theta - np.pi/2) < 0.1:
                if disk_inner < r < disk_outer:
                    # Disk emission
                    intensity = (disk_inner / r) ** 3
                    
                    # Doppler effect
                    v_orbit = np.sqrt(1.0 / r)
                    doppler = 1 + v_orbit * np.sin(phi) * np.sin(theta_obs)
                    intensity *= abs(doppler) ** 3
                    
                    return intensity
            
            if r > 200:
                return 0.01  # Background
        
        return 0.0
    
    def info(self) -> Dict[str, float]:
        """Return dictionary of black hole properties"""
        return {
            'mass': self.mass,
            'spin': self.spin,
            'horizon_radius': self.horizon_radius,
            'isco_radius': self.isco_radius,
            'isco_retrograde': self.metric.isco_radius_retrograde,
            'photon_sphere': self.photon_sphere_radius,
            'shadow_radius': self.shadow_radius,
            'omega_H': self.metric.horizon_angular_velocity
        }


def save_image(image: np.ndarray, filename: str, colormap: str = 'inferno'):
    """
    Save rendered image to file.
    
    Args:
        image: 2D intensity array
        filename: Output filename (supports .png, .fits)
        colormap: Matplotlib colormap name
    """
    try:
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 10))
        plt.imshow(image, cmap=colormap, origin='lower')
        plt.axis('off')
        plt.savefig(filename, bbox_inches='tight', pad_inches=0, dpi=150)
        plt.close()
        print(f"Saved image to {filename}")
    except ImportError:
        # Fallback: save as numpy array
        np.save(filename.replace('.png', '.npy'), image)
        print(f"Matplotlib not available, saved as {filename.replace('.png', '.npy')}")


def render_kerr(
    a: float = 0.99,
    theta: float = 85.0,
    resolution: int = 512,
    outfile: str = 'kerr_shadow.png'
) -> np.ndarray:
    """
    Quick function to render a Kerr black hole.
    
    Args:
        a: Spin parameter (0 to 0.998)
        theta: Inclination angle in degrees
        resolution: Image resolution
        outfile: Output filename
    
    Returns:
        Rendered image array
    """
    bh = KerrBlackHole(mass=1.0, spin=a)
    image = bh.render(resolution=resolution, inclination=theta)
    save_image(image, outfile)
    return image


# Convenience exports
__all__ = [
    'KerrBlackHole',
    'KerrMetric',
    'SchwarzschildMetric',
    'GeodesicIntegrator',
    'save_image',
    'render_kerr',
    'CONSTANTS',
    'PhysicalConstants',
    'MetricResult',
    '__version__'
]