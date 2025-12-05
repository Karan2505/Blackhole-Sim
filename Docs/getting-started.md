# Getting Started with BlackHole-Sim

Welcome to BlackHole-Sim! This guide will help you get up and running quickly.

## 🚀 Quick Start (No Installation Required)

The easiest way to use BlackHole-Sim is through our live demo:

**[Launch Live Demo →](https://karan2505.github.io/Blackhole-Sim/)**

## 💻 Local Installation

### Option 1: Direct Download

1. Download or clone the repository:
   ```bash
   git clone https://github.com/karan2505/Blackhole-Sim.git
   cd Blackhole-Sim
   ```

2. Open `index.html` directly in your browser, or use a local server:
   ```bash
   # Python 3
   python -m http.server 8000
   
   # Then visit http://localhost:8000
   ```

### Option 2: Docker

```bash
# Pull and run the container
docker run -p 8080:80 ghcr.io/karan2505/Blackhole-Sim:latest

# Visit http://localhost:8080
```

## 🎮 Using the Interface

### Main Controls

| Control | Function |
|---------|----------|
| **Mass (M)** | Black hole mass in geometric units |
| **Spin (a/M)** | Dimensionless spin parameter (0 = Schwarzschild, 0.998 = near-extremal Kerr) |
| **Inclination** | Observer viewing angle (0° = face-on, 90° = edge-on) |
| **Observer Distance** | Distance from black hole in units of M |

### Visualization Toggles

- **Accretion Disk** - Show/hide the accretion disk
- **Photon Sphere** - Display the photon orbit radius
- **Geodesic Rays** - Show sample light ray paths
- **Relativistic Jet** - Display bipolar jets
- **Ergosphere** - Show the ergosphere boundary (Kerr only)

### Playback Controls

- **▶/⏸** - Play/Pause animation
- **↺** - Reset camera position
- **📷** - Capture screenshot (PNG)
- **💾** - Export data (HDF5 format)
- **🎬** - Record video

## 🔬 Understanding the Physics

### Key Quantities Displayed

| Quantity | Symbol | Description |
|----------|--------|-------------|
| Schwarzschild radius | r_s | Event horizon for a=0 |
| Event horizon | r_+ | Outer horizon radius |
| ISCO radius | r_ISCO | Innermost stable circular orbit |
| Photon sphere | r_ph | Unstable photon orbit |
| Ergosphere | r_ergo | Boundary of ergosphere at equator |
| Horizon frequency | Ω_H | Angular velocity of the horizon |

### M87* Parameters

The default parameters (a/M = 0.94, inclination = 17°) are chosen to match observations of the M87* black hole by the Event Horizon Telescope.

## 📊 GRMHD Simulation

The GRMHD (General Relativistic Magnetohydrodynamics) panel shows:

1. **2D density slice** of the accretion flow
2. **Conservation diagnostics** (mass, energy, ∇·B errors)
3. **Accretion rate** in geometric units

### Running a GRMHD Simulation

1. Adjust the **Torus Inner Radius** to set initial conditions
2. Set the **Magnetic β₀** parameter
3. Click **▶ Run Simulation** to start
4. Use **⏸ Pause** to stop and analyze
5. Click **↺ Reset to IC** to restart

## 🎨 Customization

### Color Maps

Choose from several scientific colormaps:
- **Inferno** - Thermal/heat visualization
- **Plasma** - Alternative thermal
- **Viridis** - Perceptually uniform
- **AFM Hot** - Classic astronomy colormap

### Metric Types

- **Schwarzschild** - Non-rotating (a = 0)
- **Kerr** - Rotating black hole
- **Kerr-Schild** - Horizon-penetrating coordinates

## 📤 Exporting Data

### Screenshot
Click the camera icon to save a PNG image of the current view.

### HDF5 Export
Click the download icon to export simulation data including:
- Metric components
- Geodesic paths
- GRMHD field data
- Radiative transfer output

## 🐛 Troubleshooting

### Performance Issues
- Reduce browser window size
- Close other GPU-intensive applications
- Try a different browser (Chrome recommended)

### Display Issues
- Ensure WebGL is enabled in your browser
- Update graphics drivers
- Try disabling browser extensions

## 📚 Next Steps

- Read the [Physics Guide](physics.md) for theoretical background
- Explore [Example Notebooks](../examples/) for analysis scripts
- Check the [API Reference](api.md) for programmatic access

---

Questions? [Open an issue](https://github.com/karan2505/Blackhole-Sim/issues) on GitHub!
