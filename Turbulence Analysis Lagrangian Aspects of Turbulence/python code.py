
# ===================================
# Lagrangian Aspects of Turbulence
# ===================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline

data = np.load('isotropic1024_stack3.npz')
ux = data['u']  
uy = data['v']
uz = data['w']

Nx, Ny, Nz = ux.shape
L = 2 * np.pi     
dx = L / Nx
dy = L / Ny
dz = L / Nz       
nu = 0.000185     

# Extract middle plane (z=1) for each velocity component
ux_mid = ux[:, :, 1]
uy_mid = uy[:, :, 1]
uz_mid = uz[:, :, 1]
# Compute x-derivatives (periodic boundaries)
du_dx = (np.roll(ux_mid, -1, axis=0) - np.roll(ux_mid, 1, axis=0)) / (2 * dx)
dv_dx = (np.roll(uy_mid, -1, axis=0) - np.roll(uy_mid, 1, axis=0)) / (2 * dx)
dw_dx = (np.roll(uz_mid, -1, axis=0) - np.roll(uz_mid, 1, axis=0)) / (2 * dx)

# Compute y-derivatives (periodic boundaries)
du_dy = (np.roll(ux_mid, -1, axis=1) - np.roll(ux_mid, 1, axis=1)) / (2 * dy)
dv_dy = (np.roll(uy_mid, -1, axis=1) - np.roll(uy_mid, 1, axis=1)) / (2 * dy)
dw_dy = (np.roll(uz_mid, -1, axis=1) - np.roll(uz_mid, 1, axis=1)) / (2 * dy)

# Compute z-derivatives using central differencing (non-periodic)
du_dz = (ux[:, :, 2] - ux[:, :, 0]) / (2 * dz)
dv_dz = (uy[:, :, 2] - uy[:, :, 0]) / (2 * dz)
dw_dz = (uz[:, :, 2] - uz[:, :, 0]) / (2 * dz)

# using the formula for dissipation give in research paper : 
# 2000-Chacin.Cantwell-JFM-Dynamics of a low Reynolds number turbulent boundary layer

sx = (du_dx * (du_dx + du_dx) + du_dy * (du_dy + dv_dx) + du_dz * (du_dz + dw_dx))
sy = dv_dx * (dv_dx + du_dy) + dv_dy * (dv_dy + dv_dy) + dv_dz * (dv_dz + dw_dy)    
sz = dw_dx * (dw_dx + du_dz) + dw_dy * (dw_dy + dv_dz) + dw_dz * (dw_dz + dw_dz)

SS = sx + sy + sz
#print('SS.shape', np.mean(SS))
epsilon =  nu * np.mean(SS)
#print('epsilon', epsilon)


tau_eta = np.sqrt(nu / epsilon)
dt = tau_eta / 20.0
print(f"Computed τη = {tau_eta}, using Δt = {dt}")

#----------------------------------------------------------

# ans 1: Lagrangian Trajectories


Np = 20
np.random.seed(42)
x = np.random.uniform(0, L, Np)
y = np.random.uniform(0, L, Np)

T_max = 10
num_steps = int(T_max / dt)



def bilinear_interp(u, x_p, y_p):
    """Bilinear interpolation for off-grid velocities."""
    x_p_mod = x_p % L
    y_p_mod = y_p % L
    i = x_p_mod * Nx / L
    j = y_p_mod * Nx / L
    i0, j0 = int(i) % Nx, int(j) % Nx
    i1, j1 = (i0 + 1) % Nx, (j0 + 1) % Nx
    di, dj = i - i0, j - j0
    return (u[i0, j0] * (1 - di) * (1 - dj) +
            u[i1, j0] * di * (1 - dj) +
            u[i0, j1] * (1 - di) * dj +
            u[i1, j1] * di * dj)


def simulate_trajectories(T_max, dt, x, y, ux_mid, uy_mid):
    num_steps = int(T_max / dt)
    x_traj = np.zeros((num_steps + 1, len(x)))
    y_traj = np.zeros((num_steps + 1, len(y)))
    x_traj[0], y_traj[0] = x, y
    for step in range(num_steps):
        for p in range(len(x)):
            u = bilinear_interp(ux_mid, x_traj[step, p], y_traj[step, p])
            v = bilinear_interp(uy_mid, x_traj[step, p], y_traj[step, p])
           
            x_traj[step+1, p] = x_traj[step, p] + u * dt
            y_traj[step+1, p] = y_traj[step, p] + v * dt
    return x_traj, y_traj

# Plot for T = 1, 5, 10
for T in [1, 5, 10]:
    x_traj, y_traj = simulate_trajectories(T, dt, x, y, ux_mid, uy_mid)
    plt.figure(figsize=(8, 6))
    for p in range(Np):
        plt.plot(x_traj[:, p], y_traj[:, p], linewidth=0.5)
    plt.title(f'Trajectories (T={T})')
    plt.xlabel('x'); plt.ylabel('y')
    plt.xlim(0, L); plt.ylim(0, L)
    plt.show()

#########################################################################
# ans 3 -----------------------------------------
#Plot the mean-square-displacement as a function of (pseudo)time

Np_msd = 10000
x_msd = np.random.uniform(0, L, Np_msd)
y_msd = np.random.uniform(0, L, Np_msd)
msd = np.zeros(num_steps + 1)
x0_msd, y0_msd = x_msd.copy(), y_msd.copy()

for step in range(num_steps):
    for p in range(Np_msd):
        u = bilinear_interp(ux_mid, x_msd[p], y_msd[p])
        v = bilinear_interp(uy_mid, x_msd[p], y_msd[p])
        x_msd[p] += u * dt
        y_msd[p] += v * dt
    msd[step+1] = np.mean((x_msd - x0_msd)**2 + (y_msd - y0_msd)**2)

plt.loglog(np.arange(num_steps+1)*dt, msd)
plt.xlabel('Time'); plt.ylabel('MSD')
plt.title('Mean-Square Displacement')
plt.show()

#########################################################################
# ans  4: Turbulent Diffusivity

t_diffusive = np.arange(num_steps+1)*dt > 5 * tau_eta
D_T = np.polyfit(np.arange(num_steps+1)[t_diffusive]*dt, msd[t_diffusive], 1)[0]
print(f'Turbulent Diffusivity: {D_T} ')  #(vs. molecular ~1e-9 m²/s)


# Time array normalized by τη
t_normalized = np.arange(num_steps + 1) * dt / tau_eta

plt.figure(figsize=(10, 6))
plt.loglog(t_normalized, msd, 'y-', linewidth=2, label='MSD')

# Add theoretical scaling lines
# 1. Ballistic regime (t^2)
t_ballistic = np.linspace(0.01, 0.1 * tau_eta, 10) / tau_eta  # Early times
C_ballistic = msd[1] / (t_normalized[1]**2)  # Match initial slope
plt.loglog(t_ballistic, C_ballistic * t_ballistic**2, 'r--', 
           label=r'Ballistic ($\sim t^2$)')

# 2. Diffusive regime (t)
t_diffusive = np.linspace(10 * tau_eta, max(t_normalized), 10) / tau_eta  # Late times
C_diffusive = D_T * tau_eta  # From previous calculation
plt.loglog(t_diffusive, C_diffusive * t_diffusive, 'g--', 
           label=r'Diffusive ($\sim t$)')

# Mark transition regions
plt.axvline(x=1, color='k', linestyle=':', alpha=0.5, label=r'$t \approx \tau_\eta$')
plt.text(0.3, 0.01, 'Ballistic\nRegime', ha='center', va='center', 
         fontsize=10, bbox=dict(facecolor='white', alpha=0.8))
plt.text(10, 0.1, 'Diffusive\nRegime', ha='center', va='center', 
         fontsize=10, bbox=dict(facecolor='white', alpha=0.8))

plt.xlabel(r'$t/\tau_\eta$', fontsize=12)
plt.ylabel(r'$\langle \Delta r^2(t) \rangle$', fontsize=12)
plt.title('MSD with Ballistic-to-Diffusive Transition', fontsize=14)
plt.legend(fontsize=10)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout()
plt.show()




#########################################################################
# ans  5: Richardson Pair Dispersion



def integrate_particles(x0, y0, n_steps, dt):
    """
    Integrate trajectories for a given initial condition using bilinear interpolation.
    x0, y0: initial positions (arrays)
    n_steps: number of time steps
    dt: time step size
    Returns: trajectories x(t), y(t) of shape [n_steps, number of particles]
    """
    Np = len(x0)
    x_traj = np.zeros((n_steps, Np))
    y_traj = np.zeros((n_steps, Np))
    x_traj[0] = x0
    y_traj[0] = y0

    for t in range(1, n_steps):
        for p in range(Np):
            x_mod = x_traj[t-1, p] % L
            y_mod = y_traj[t-1, p] % L
            u_part = bilinear_interp(ux_mid, x_mod, y_mod)
            v_part = bilinear_interp(uy_mid, x_mod, y_mod)

            x_traj[t, p] = x_traj[t-1, p] + u_part * dt
            y_traj[t, p] = y_traj[t-1, p] + v_part * dt

    return x_traj, y_traj

def integrate_pair_trajectories(xA0, yA0, eps, n_steps, dt):
    """
    Initialize Set B by perturbing Set A with fixed separation eps in random directions.
    Returns trajectories for both sets: (xA, yA) and (xB, yB).
    """
    Np = len(xA0)
    theta = np.random.uniform(0, 2*np.pi, size=Np)
    xB0 = xA0 + eps * np.cos(theta)
    yB0 = yA0 + eps * np.sin(theta)
    
    xA, yA = integrate_particles(xA0, yA0, n_steps, dt)
    xB, yB = integrate_particles(xB0, yB0, n_steps, dt)
    return xA, yA, xB, yB

#########################################################################
# (a) Plot pair trajectories for Np = 20, ε = 0.5Δx, T = 10.

Np_pair = 20
eps = 0.5 * dx
particles_xA0 = np.random.uniform(0, L, size=Np_pair)
particles_yA0 = np.random.uniform(0, L, size=Np_pair)
n_steps_pair = int(10 / dt)
xA, yA, xB, yB = integrate_pair_trajectories(particles_xA0, particles_yA0, eps, n_steps_pair, dt)


plt.figure(figsize=(6, 6))

# Plot first point of each set with label
plt.plot(xA[0], yA[0], '.', color='k', label='Set A')
plt.plot(xB[0], yB[0], '.', color='r', label='Set B')

# Plot the rest without label
plt.plot(xA, yA, '.', linestyle='None', color='k')
plt.plot(xB, yB, '.', linestyle='None', color='r')

plt.xlabel('x')
plt.ylabel('y')
plt.title('Pair Trajectories for Np = 20, ε = 0.5Δx, T = 10')
plt.legend()
plt.xlim(0, L)
plt.ylim(0, L)
plt.grid(True)
plt.gca().set_aspect('equal')
plt.show()

#########################################################################
# (b) Log-log plot of pair separation (⟨Δr²(t)⟩/ε²) vs t/τ_η for various ε.

eps_list = [0.1*dx, 0.5*dx, 1.0*dx, 5.0*dx, 10.0*dx]
Np_dispersion = 100  # Adjust if needed; ideally up to 10^4 if possible
T_disp = 10
n_steps_disp = int(T_disp / dt)
time_disp = np.arange(n_steps_disp) * dt

plt.figure()
for eps_val in eps_list:
    # Initialize Set A for pair dispersion
    xA0 = np.random.uniform(0, L, size=Np_dispersion)
    yA0 = np.random.uniform(0, L, size=Np_dispersion)
    # Integrate pair trajectories
    xA, yA, xB, yB = integrate_pair_trajectories(xA0, yA0, eps_val, n_steps_disp, dt)
    # Compute separation squared for each pair at each time
    separation_sq = (xA - xB)**2 + (yA - yB)**2
    mean_sep_sq = np.mean(separation_sq, axis=1)
    plt.loglog(time_disp / tau_eta, mean_sep_sq / (eps_val**2), label=f'ε = {eps_val:.2e}')
    
    # Optionally, perform a line fit in the inertial range (here, using a subset of time indices)
    # Uncomment and adjust indices if needed:
    # inertial_idx = (time_disp/tau_eta > 1) & (time_disp/tau_eta < 5)
    # p_fit = np.polyfit(np.log(time_disp[inertial_idx]), np.log(mean_sep_sq[inertial_idx]), 1)
    # plt.loglog(time_disp[inertial_idx] / tau_eta, np.exp(np.polyval(p_fit, np.log(time_disp[inertial_idx])))/ (eps_val**2), '--')
    
plt.xlabel('t/τ_η')
plt.ylabel('⟨Δr²(t)⟩/ε²')
plt.title('Richardson Pair Dispersion (Log-Log)')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.show()



# Pair Separation Curves for Specific ε/Δx Values

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define power-law function for fitting
def power_law(t, a, b):
    return a * t**b

# Parameters
eps_list2 = [0.1, 0.5, 1.0, 5.0, 10]  # ε/Δx values (as specified)
Np_dispersion = 100                # Number of particle pairs
dx = L / Nx                            # Grid spacing (Δx)

# Compute Kolmogorov and integral scales
#epsilon = 2 * nu * np.mean(S_squared)  # Dissipation rate
tau_eta = np.sqrt(nu / epsilon)        # Kolmogorov time
tau_L = 1.364 / np.sqrt(np.mean(ux_mid**2 + uy_mid**2 + uz_mid**2))  # Integral time

# Inertial range mask (t/τη ≈ 1 to τ_L/τ_η)
inertial_mask = (time_disp >= tau_eta) & (time_disp <=  tau_L)

# Plot settings
colors = ['b', 'g', 'r', 'c', 'm']    # Colors for each ε

# Individual plots for each ε/Δx
for eps_val, color in zip(eps_list2, colors):
    plt.figure(figsize=(10, 6))
    
    # Initialize particle pairs with ε = eps_val * Δx
    xA0 = np.random.uniform(0, L, size=Np_dispersion)
    yA0 = np.random.uniform(0, L, size=Np_dispersion)
    theta = np.random.uniform(0, 2*np.pi, size=Np_dispersion)
    xB0 = xA0 + (eps_val * dx) * np.cos(theta)
    yB0 = yA0 + (eps_val * dx) * np.sin(theta)
    
    # Simulate trajectories (use your existing function)
    xA, yA, xB, yB = integrate_pair_trajectories(xA0, yA0, eps_val * dx, n_steps_disp, dt)
    separation_sq = (xA - xB)**2 + (yA - yB)**2
    mean_sep_sq = np.mean(separation_sq, axis=1)
    
    # Normalize and mask inertial range
    t_inertial = time_disp[inertial_mask] / tau_eta
    sep_inertial = mean_sep_sq[inertial_mask] / (eps_val * dx)**2
    
    # Plot data points
    plt.loglog(t_inertial, sep_inertial, 'o', color=color, markersize=4, 
               label=f'ε = {eps_val}Δx')
    
    # Fit power law (Δr² ∼ t^b)
    params, _ = curve_fit(power_law, t_inertial, sep_inertial, p0=[0.1, 3])
    fitted_sep = power_law(t_inertial, *params)
    
    # Plot fitted curve
    plt.loglog(t_inertial, fitted_sep, 'k-', linewidth=1.5,
              label=rf'Fit: $t^{{{params[1]:.2f}}}$')
    
    # Theoretical t³ line (adjust prefactor to match data)
    plt.loglog(t_inertial, 0.05 * t_inertial**3, 'r--', 
               label=r'Theory: $t^3$', linewidth=1.5)
    
    # Formatting
    plt.xlabel(r'$t / \tau_\eta$', fontsize=12)
    plt.ylabel(r'$\langle \Delta r^2 \rangle / \epsilon^2$', fontsize=12)
    plt.title(f'Pair Separation (ε = {eps_val}Δx)', fontsize=14)
    plt.legend(fontsize=10, loc='upper left')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.xlim(1, 10)
    plt.ylim(1e-2, 1e4)
    plt.tight_layout()
    plt.show()
    
    print(f"ε = {eps_val}Δx | Fitted exponent: {params[1]:.2f} (expected: 3)")

#########################################################################
# (c) Semilog plot of pair separation to observe early exponential growth.

plt.figure()
for eps_val in eps_list:
    # Initialize new pairs
    xA0 = np.random.uniform(0, L, size=Np_dispersion)
    yA0 = np.random.uniform(0, L, size=Np_dispersion)
    xA, yA, xB, yB = integrate_pair_trajectories(xA0, yA0, eps_val, n_steps_disp, dt)
    separation_sq = (xA - xB)**2 + (yA - yB)**2
    mean_sep_sq = np.mean(separation_sq, axis=1)
    
    # Normalize time by Kolmogorov time scale
    time_normalized = time_disp / tau_eta
    plt.semilogy(time_normalized, mean_sep_sq / (eps_val**2), 
                label=f'ε = {eps_val:.2e}')

plt.xlabel('t/τ_η')  # Changed from 'Time' to normalized time
plt.ylabel('⟨Δr²(t)⟩/ε² (log scale)')
plt.title('Pair Separation (Semilog) - Early Exponential Growth')
plt.legend()
plt.grid(True)
plt.show()

# Lyapunov exponent calculation with normalized time
eps_fit = eps_list[1]  # choose, say, ε = 0.5Δx
xA0 = np.random.uniform(0, L, size=Np_dispersion)
yA0 = np.random.uniform(0, L, size=Np_dispersion)
xA, yA, xB, yB = integrate_pair_trajectories(xA0, yA0, eps_fit, n_steps_disp, dt)
separation_sq = (xA - xB)**2 + (yA - yB)**2
mean_sep_sq = np.mean(separation_sq, axis=1)

# Normalize time for fitting
time_normalized = time_disp / tau_eta

# Select early time window for fitting (first 10% of simulation)
fit_end = int(0.1 * n_steps_disp)
t_fit = time_normalized[:fit_end]
log_sep = np.log(mean_sep_sq[:fit_end])

# Fit line: log(⟨Δr²⟩) = 2Λ*(t/τ_η) + constant
# Note the factor of 2 because Δr² ∝ exp(2Λt)
fit_coef = np.polyfit(t_fit, log_sep, 1)
Lambda = fit_coef[0]  # Extract Λ from the slope
print(f"Estimated Lyapunov exponent Λ ≈ {Lambda:.4f}")

# Plot results with normalized time axis
plt.figure()
plt.semilogy(time_normalized, mean_sep_sq / (eps_fit**2), 'b.', label='Data')
plt.semilogy(t_fit, np.exp(np.polyval(fit_coef, t_fit)) / (eps_fit**2), 'r-', lw=2,
             label=f'Fit: Λ = {Lambda:.4f}')
plt.xlabel('t/τ_η')  # Changed to normalized time
plt.ylabel('⟨Δr²(t)⟩/ε² (log scale)')
plt.title('Early Time Exponential Growth (Normalized Time)')
plt.legend()
plt.grid(True)
plt.show()