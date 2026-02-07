# Lid-Driven Cavity Flow using SIMPLE (Staggered Grid)

This repository contains a 2D finite-volume solver for the **incompressible, laminar, viscous Navier–Stokes equations**  
using the **SIMPLE algorithm** on a **staggered grid**, applied to the classical lid-driven cavity problem.

---

## Problem Description
- Flow type: Incompressible, laminar, constant density
- Geometry: Square cavity
- Top lid velocity: 1.0
- All other walls: No-slip
- Pressure–velocity coupling: SIMPLE method
- Discretization: Finite Volume Method (central differencing)

---

## Recommended Grid Sizes
| Reynolds Number | Grid Size |
|-----------------|-----------|
| Re = 100        | 51 × 51   |
| Re = 1000       | 501 × 501 |

> Note: Higher Reynolds numbers require finer grids for numerical stability and accuracy.

---

## Parameters Used
- Kinematic viscosity (`nu_f`): 0.01 → 0.001  
- Density (`rho_f`): 1.0  
- Lid velocity (`U`): 1.0  
- Domain length (`L`): 1.0  

You may change these parameters to test different Reynolds numbers and verify the solution behavior.

---

## Post-processing & Visualization
Python scripts are provided to visualize:
- Velocity components (u, v)
- Pressure field
- Centerline velocity profiles

The visualization code was generated with assistance from ChatGPT and is intended for post-processing only.

You can run the Python scripts using Jupyter Notebook:
https://jupyter.org/try-jupyter/lab/

---

## Notes
- This solver is intended for **educational and research purposes**
- The code has not been optimized for performance
- Stability at high Reynolds numbers may require finer grids or improved discretization schemes

---

## Issues & Contributions
Feel free to **raise an issue** if you encounter any problems or have suggestions for improvement.
