# Contributing to BlackHole-Sim

First off, thank you for considering contributing to BlackHole-Sim! 🌌

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [How to Contribute](#how-to-contribute)
- [Development Setup](#development-setup)
- [Pull Request Process](#pull-request-process)
- [Style Guidelines](#style-guidelines)

## Code of Conduct

This project and everyone participating in it is governed by our commitment to providing a welcoming and inclusive environment. Please be respectful and constructive in all interactions.

## Getting Started

### Prerequisites

- Git
- A modern web browser (Chrome, Firefox, Safari, Edge)
- Optional: Python 3.8+ for local development server
- Optional: Node.js for advanced development

### Quick Setup

```bash
# Fork and clone the repository
git clone https://github.com/karan2505/Blackhole-Sim.git
cd Blackhole-Sim

# Start a local server
python -m http.server 8000

# Visit http://localhost:8000 in your browser
```

## How to Contribute

### Reporting Bugs

Before creating bug reports, please check existing issues. When creating a bug report, include:

- **Clear title** describing the issue
- **Steps to reproduce** the behavior
- **Expected behavior** vs actual behavior
- **Screenshots** if applicable
- **Browser and OS** information

### Suggesting Features

Feature requests are welcome! Please provide:

- **Clear description** of the feature
- **Use case** - why would this be useful?
- **Possible implementation** approach (optional)

### Code Contributions

1. **Find an issue** to work on or create one
2. **Comment** on the issue to let others know you're working on it
3. **Fork** the repository
4. **Create a branch** for your feature
5. **Make your changes**
6. **Test** your changes locally
7. **Submit** a pull request

## Development Setup

### Project Structure

```
Blackhole-Sim/
├── index.html          # Main application (single-page)
├── README.md           # Project documentation
├── LICENSE             # MIT License
├── CONTRIBUTING.md     # This file
├── .github/
│   └── workflows/
│       └── deploy.yml  # CI/CD pipeline
└── docs/               # Additional documentation
```

### Local Development

```bash
# Using Python's built-in server
python -m http.server 8000

# Using Node.js (if you have it)
npx serve .

# Using PHP (if you have it)
php -S localhost:8000
```

### Testing Changes

1. Open the browser developer tools (F12)
2. Check the Console for any JavaScript errors
3. Test all interactive features:
   - Parameter sliders
   - Visualization toggles
   - Play/pause controls
   - Export functions

## Pull Request Process

1. **Update documentation** if you change functionality
2. **Test thoroughly** on multiple browsers if possible
3. **Write a clear PR description** explaining:
   - What changes you made
   - Why you made them
   - Any testing you performed
4. **Link related issues** using keywords like "Fixes #123"
5. **Request review** from maintainers

### PR Title Format

Use conventional commit style:
- `feat: add new visualization mode`
- `fix: correct geodesic calculation`
- `docs: update physics equations`
- `style: improve mobile responsiveness`
- `refactor: optimize render loop`

## Style Guidelines

### HTML/CSS

- Use **Tailwind CSS** utility classes when possible
- Keep custom CSS minimal and well-organized
- Use semantic HTML elements
- Maintain accessibility (ARIA labels, alt text)

### JavaScript

- Use **ES6+** features (const/let, arrow functions, etc.)
- Add comments for complex physics calculations
- Keep functions focused and modular
- Handle errors gracefully

### Example Code Style

```javascript
// Good: Clear, documented, modular
/**
 * Calculate the ISCO radius for a Kerr black hole
 * @param {number} M - Black hole mass
 * @param {number} a - Spin parameter (|a| < M)
 * @returns {number} ISCO radius in units of M
 */
function calculateISCO(M, a) {
    const z1 = 1 + Math.pow(1 - a*a/(M*M), 1/3) * 
               (Math.pow(1 + a/M, 1/3) + Math.pow(1 - a/M, 1/3));
    const z2 = Math.sqrt(3*a*a/(M*M) + z1*z1);
    return M * (3 + z2 - Math.sqrt((3 - z1) * (3 + z1 + 2*z2)));
}
```

### Physics Accuracy

When modifying physics calculations:
- **Cite sources** for equations
- **Include units** in comments
- **Validate** against known analytic solutions
- **Document limitations** and approximations

## Questions?

Feel free to:
- Open an issue with the `question` label
- Start a discussion in the Discussions tab
- Reach out to maintainers

---

Thank you for contributing to BlackHole-Sim! 🚀
