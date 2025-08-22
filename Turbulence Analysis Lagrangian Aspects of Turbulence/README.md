
---

### **Lagrangian Aspects of Turbulence** 🧭

This section explores turbulence from a Lagrangian perspective by tracking the motion of tracer particles in a frozen velocity field using Euler integration.

#### 🔹 Lagrangian Trajectories

* Particle paths are plotted for different counts: small (**Np = 20**) and large (**Np = 1000**).
* For long simulation times, particles are seen to **circulate, drift, or become trapped** in vortex structures.
* Dense regions in the velocity field lead to recirculation and trapping of some particles.

#### 🔹 Mean-Square Displacement (MSD)

* MSD is plotted on a **log-log scale** to identify different dispersion regimes:

  * **Ballistic regime** at short times (MSD ∝ t²).
  * **Diffusive regime** at long times (MSD ∝ t).
* The computed **turbulent diffusivity** is significantly greater than molecular diffusivity, showing how turbulence enhances mixing.

#### 🔹 Richardson Pair Dispersion

* The separation of initially close particle pairs is analyzed over time.
* Trajectories of such pairs diverge due to chaotic flow behavior.

#### 🔹 Scaling Laws & Lyapunov Exponent

* **Normalized pair separation vs. time** is plotted on both **log-log** and **semi-log** scales.
* Log-log plots verify **Richardson’s t³ scaling** in the inertial range.
* Semi-log plots reveal **early exponential separation**, from which the **Lyapunov exponent** is calculated — quantifying the chaos level in particle dispersion.

---

