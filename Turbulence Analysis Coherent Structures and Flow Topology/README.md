
---

## Turbulence Analysis: Coherent Structures and Lagrangian Dynamics 🌪️🧭

This report, based on a computational assignment, presents a comprehensive analysis of two fundamental aspects of turbulence — **coherent structures** and **Lagrangian particle dynamics** — using a dataset from a Direct Numerical Simulation (DNS) of homogeneous, isotropic turbulence.

---

### **Part 1: Coherent Structures and Flow Topology** 🌀

This section investigates the local structure of turbulent flows by analyzing the **velocity gradient tensor**.

#### 🔹 Invariants and the QR Plane

* The velocity gradient tensor is characterized by three invariants: **P**, **Q**, and **R**.
* Due to fluid incompressibility, **P is approximately zero**, and the Q-R plane becomes the focus for classifying local flow topologies.
* The characteristic equation is derived, and the relation between the invariants and strain-rate and rotation-rate tensors is discussed.

#### 🔹 Eigenvalue Analysis

* Eigenvalues of the velocity gradient tensor are computed for each grid point in a 2D plane.
* The probability density functions (PDFs) of these eigenvalues are plotted on a semi-log scale, revealing statistical trends in the flow dynamics.

#### 🔹 Incompressibility Verification

* The mean of the absolute value of the **P-invariant** is confirmed to be near zero, validating the incompressibility assumption.

#### 🔹 Q and R Fields

* Normalized spatial fields of **Q** and **R** are plotted.
* Both fields exhibit patchy, filament-like structures, typical in turbulence.
* High Q regions correlate with rotational-dominant flow, and the sign of R indicates **vortex stretching** (R > 0) or **vortex compression** (R < 0).

#### 🔹 Joint Distribution (Q vs. R)

* A joint probability distribution of Q and R is plotted.
* The result shows the classic **“tear-drop” shape**, a universal feature in turbulent flows.
* Dominant flow topologies are identified: **vortex stretching** and **biaxial strain**.
* The analysis links these structures to the **enstrophy equation**, showing how vortex stretching enhances enstrophy, while biaxial strain reduces it.

---

