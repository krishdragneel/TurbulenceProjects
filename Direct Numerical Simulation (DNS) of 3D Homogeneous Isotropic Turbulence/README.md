
---

# 📘 Direct Numerical Simulation (DNS) of 3D Homogeneous Isotropic Turbulence

This report presents a detailed study of **Direct Numerical Simulation (DNS)** of **three-dimensional homogeneous, isotropic turbulence**. The work includes both **theoretical derivations** and **computational analysis** of a 2D slice of a high-resolution turbulence dataset, aiming to understand the dynamics and structure of turbulent flows.

---

## 🧪 Part 1: Fluid Flow Equations

This section focuses on deriving and interpreting key equations from the incompressible Navier-Stokes equations.

### 🔹 Kinetic Energy Equation

* The transport equation for **kinetic energy** is derived by manipulating the Navier-Stokes equations.
* In regions far from walls (unbounded turbulence), several terms become negligible—such as convective and viscous flux surface integrals, and pressure work.
* The simplified form of the equation reveals that **kinetic energy is dissipated by viscosity**.

### 🔹 Enstrophy Equation

* The **enstrophy** transport equation is derived using the curl of the Navier-Stokes equations, followed by a projection onto the vorticity field.
* While kinetic energy is redistributed via pressure, **enstrophy grows due to vortex stretching** and dissipates through vorticity gradients.
* Unlike kinetic energy, **pressure does not directly affect enstrophy**.

---

## 💻 Part 2: Turbulence Data Analysis

This part analyzes a **1024 × 1024** grid turbulence dataset using Python and numerical tools.

### 🔹 Flow Characteristics

* The **Reynolds number** is calculated to be approximately **4478**, indicating a highly turbulent regime.
* The **Kolmogorov scales** are computed, and the grid resolution is found to be **about 2.5 times coarser** than the Kolmogorov length scale.

### 🔹 Field Visualizations

* **Kinetic energy fields** display smooth, filament-like structures.
* **Vorticity fields** show complex, small-scale, and alternating patterns—capturing the intermittent nature of turbulence.

### 🔹 Statistical Analysis

* **Probability Density Functions (PDFs)** and **Cumulative Distribution Functions (CDFs)** are generated for:

  * Velocity components
  * Velocity gradients
  * Vorticity
  * Enstrophy
* Velocity component PDFs follow a **Gaussian distribution**.
* PDFs for velocity gradients and vorticity **deviate significantly from Gaussian**, indicating skewness, high kurtosis, and **intermittency**—a hallmark of turbulence.
* Enstrophy PDFs are **heavily skewed positive**, consistent with it being a squared quantity.

### 🔹 Local Reynolds Number Field

* A spatial map of the **local Reynolds number** is computed.
* Regions with high kinetic energy **often (but not always)** align with regions of high local Reynolds number.
* **Vorticity** strongly influences the localized Reynolds number, especially around small-scale eddies and vortices.
* The local Reynolds number field appears **patchy and less smooth** than the energy and vorticity fields, highlighting fine-scale variation.

---

