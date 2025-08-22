
---

## 💨 **Turbulence Projects Summary**

---

### **Project 1: Fluid Flow Equations & Data Analysis** 🧪

This project combined **theoretical derivations** and **computational data analysis** to investigate turbulent flows.

#### 🔹 Theoretical Derivations

* Derived the **kinetic energy** and **enstrophy** transport equations from the incompressible Navier-Stokes equations.
* **Kinetic energy** is dissipated through **velocity gradients**, while **enstrophy** is dissipated through **vorticity gradients**.
* **Enstrophy** production occurs via **vortex stretching**, whereas **kinetic energy** is redistributed by **pressure effects**.

#### 🔹 Data Analysis

* Used a **1024×1024 DNS dataset** to compute the **Reynolds number** and **Kolmogorov scales**.
* Plotted normalized **kinetic energy** and **vorticity fields**, revealing distinct turbulent structures.
* Conducted **statistical analysis**:

  * Probability Density Functions (PDFs) and Cumulative Distribution Functions (CDFs) were plotted for velocity components, velocity gradients, and vorticity.
  * Results deviated from Gaussian distributions, indicating **intermittency** in turbulence.
* Plotted a **local Reynolds number field**, which correlated with but did not exactly match the kinetic energy and vorticity fields.

---

### **Project 2: Turbulence Scaling Laws** 📈

This project focused on verifying **turbulence scaling theories** using the same **1024×1024 dataset**.

#### 🔹 Fourier Space Analysis

* Verified **Parseval’s theorem**, confirming conservation of total energy between real and Fourier space.
* Plotted the **1D and 2D energy spectra**:

  * The **1D spectrum** followed **Kolmogorov’s -5/3 scaling law** in the inertial range.
  * The **2D shell-averaged spectrum** also showed scaling behavior, though with slight deviation in the exponent.

#### 🔹 Real Space Analysis

* Plotted **longitudinal and transverse velocity correlation functions**, showing how velocity correlations decay with distance.
* Analyzed **longitudinal structure functions** for different orders:

  * **Kolmogorov’s 4/5ths law** was confirmed for the third-order structure function.
* Used **Extended Self-Similarity (ESS)** plots to extract scaling exponents (ζₚ).

  * Found deviations from the theoretical linear scaling (p/3), confirming **anomalous scaling** due to **intermittent energy dissipation**.

---

### **Project 3: Coherent Structures & Lagrangian Dynamics** 🌪️🧭

This project analyzed **coherent flow structures** and **Lagrangian particle motion** using a **3D turbulence dataset**.

#### 🔹 Flow Topology and Coherent Structures

* Investigated the **velocity gradient tensor** and its **invariants** (P, Q, R).

  * Verified **incompressibility** by confirming that the P-invariant was near zero.
  * Described flow topology using the **Q-R plane**, which showed the classic **“tear-drop” shape**.
* Identified dominant structures such as:

  * **Vortex stretching** (R > 0, Q > 0)
  * **Biaxial strain** (R < 0, Q < 0)
* Plotted **Q and R fields**, which revealed **filamentary and patchy flow structures** typical of turbulence.

#### 🔹 Lagrangian Particle Dynamics

* Simulated **tracer particle trajectories** in a **frozen velocity field**:

  * Observed particles becoming **trapped in vortex structures** or recirculation zones.
* Plotted **Mean-Square Displacement (MSD)**:

  * Identified a transition from **ballistic motion** (\~t²) at short times to **diffusive behavior** (\~t) at longer times.
* Analyzed **Richardson pair dispersion**:

  * Tracked how particle pairs initially close together separated over time.
  * Detected an **exponential growth regime**, from which the **Lyapunov exponent** was extracted — a measure of **chaotic separation**.

---

