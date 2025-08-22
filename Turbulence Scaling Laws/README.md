
---

# 🌪️ Turbulence Scaling Laws: Analysis from DNS Data

This project presents a comprehensive analysis of turbulence scaling laws using a **2D dataset from a Direct Numerical Simulation (DNS)** of homogeneous, isotropic turbulence. The study validates core theoretical principles in both **Fourier space** and **real space**, with a focus on energy transfer across scales and the statistical nature of turbulence.

---

## 📊 Part 1: Fourier Space Analysis

### 🔹 Parseval's Theorem

* The project verifies **Parseval’s theorem**, confirming that total energy is conserved between real and Fourier domains.
* Total energy in:

  * **Real space**: 428.7784573503308
  * **Fourier space**: 428.7784573503309
* The near-identical values confirm the numerical accuracy of the Fourier transform used in spectral analysis.

### 🔹 Energy Spectrum

* The **1D energy spectrum** (averaged over the y-direction) is plotted on a log-log scale.
* The plot shows a **power-law scaling** consistent with **Kolmogorov's -5/3 law** in the inertial range, reflecting the expected energy cascade in turbulence.
* The **2D energy spectrum** is also computed using **shell averaging** (averaging over constant wavenumber shells).

  * The fitted exponent from this method is approximately **-2.68**, deviating from the theoretical -5/3 due to dataset resolution or numerical artifacts.

---

## 🌀 Part 2: Real Space Analysis

### 🔹 Velocity Correlation Functions

* **Longitudinal** and **transverse velocity correlations** are plotted to examine the spatial coherence of velocity fluctuations.
* Strong correlation is observed at short distances, decreasing with separation.
* Negative correlations at certain distances imply flow structures with **opposing velocity directions**, a signature of turbulent eddies.

### 🔹 Longitudinal Structure Functions

* The **longitudinal structure functions**, which measure the velocity difference between points separated by distance *r*, are plotted for orders **p = 1 to 7** on a log-log scale.

### 🔹 Kolmogorov’s 4/5 Law

* The **third-order structure function** is examined to verify **Kolmogorov’s 4/5 law**, a key result in turbulence theory.
* The analysis confirms the expected linear scaling behavior in the inertial range.

### 🔹 Extended Self-Similarity (ESS)

* To enhance the scaling range, the **ESS technique** is used, where higher-order structure functions are plotted against the third-order function.
* This method provides more robust scaling behavior even when the inertial range is limited.

### 🔹 Scaling Exponents

* **Scaling exponents** for each structure function order are computed and compared to Kolmogorov's prediction (p/3).
* Significant deviations from this prediction are observed, especially at higher orders.
* This **anomalous scaling** is attributed to **intermittency**—irregular, intense fluctuations in energy dissipation at small scales.
* Higher-order structure functions are more sensitive to these intermittent events, explaining the growing deviation with increasing p.

---

✅ **Conclusion**
The project successfully validates key aspects of turbulence theory:

* Conservation of energy across domains (Parseval's theorem)
* Scale-dependent energy distribution (Kolmogorov scaling)
* Statistical structure of velocity differences and intermittency effects
  It demonstrates both the **power and limits** of theoretical models when applied to real DNS data.

---


