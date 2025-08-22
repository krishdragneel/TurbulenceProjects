import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate
from scipy.stats import linregress

# Load dataset
file_path = "isotropic1024_slice.npz"
data = np.load(file_path)

# Extract velocity fields
u = data['u']  # Velocity in x-direction
v = data['v']  # Velocity in y-direction

###############################
# 1. Verify Parseval's Theorem #
###############################

y_index = 0  # Line cut at y=0
u_x = u[y_index, :]

# Compute Fourier Transform
u_x_hat = np.fft.fft(u_x)
power_spectral_density = np.abs(u_x_hat) ** 2

# Compute energy in real and Fourier space
energy_real = np.sum(u_x**2)
energy_fourier = np.sum(power_spectral_density) / len(u_x)

print("Energy in real space:", energy_real)
print("Energy in Fourier space:", energy_fourier)

u_hat = np.fft.fft(u, axis=1)

# Energy spectrum calculation
E_k_x = (np.abs(u_hat) ** 2 ) / 2
E_k_x_avg = np.mean(E_k_x, axis=0)

# Compute wavenumbers and distances
Nx = u.shape[1]
Lx = 2 * np.pi
k_x = np.fft.fftfreq(Nx, d=Lx / Nx) * (2 * np.pi)
k_x = np.fft.fftshift(k_x)
E_k_x_avg = np.fft.fftshift(E_k_x_avg)
r_values = np.arange(Nx) * (Lx / Nx)

# Plot energy spectra side by side
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Real space energy spectrum
axs[0].loglog(r_values, u_x**2, label="Real Space", color='red')
axs[0].set_xlabel("Distance r")
axs[0].set_ylabel("Energy Spectrum E(r)")
axs[0].legend()
axs[0].set_title("Real Space Energy Spectrum")

# Fourier space energy spectrum
axs[1].loglog(k_x[Nx//2+1:], E_k_x_avg[Nx//2+1:], label="Fourier Space", color='red')
axs[1].set_xlabel("Wavenumber k")
axs[1].set_ylabel("Energy Spectrum E(k)")
axs[1].legend()
axs[1].set_title("Fourier Space Energy Spectrum")

plt.tight_layout()
plt.show()

#####################################
# 2. Compute 1D Energy Spectrum E(k) #
#####################################

# Compute Fourier Transform along x for all y values
u_hat = np.fft.fft(u, axis=1)
v_hat = np.fft.fft(v, axis=1)

# Energy spectrum calculation
E_k_x = (np.abs(u_hat) ** 2 + np.abs(v_hat) ** 2) / 2
E_k_x_avg = np.mean(E_k_x, axis=0)

# Print the averaged energy spectrum values
print("Total Energy Spectrum Along x-direction Averaged Over y:")
print(np.sum(E_k_x_avg))

# Parameters
L = 1.364
nu = 0.000185
u_rms = np.sqrt(np.mean(u**2 + v**2))/np.sqrt(2)

epsilon = (u_rms**3) / L
eta = (nu**3 / epsilon)**(1/4)

# print('u_rms', u_rms)
# print('epsilon', epsilon)
# print('eta', eta)


# Inertial range calculation
inertial_range_min = 60 * eta
inertial_range_max = L / 6

# Convert to wavenumbers
k_inertial_min = 2 * np.pi / inertial_range_max
k_inertial_max = 2 * np.pi / inertial_range_min

# print(f"Inertial range: {inertial_range_min:.6f} to {inertial_range_max:.6f}")
# print(f"Wavenumber range: {k_inertial_min:.6f} to {k_inertial_max:.6f}")

# Compute wavenumbers
Nx = u.shape[1]
Lx = 2 * np.pi
k_x = np.fft.fftfreq(Nx, d=Lx / Nx) * (2 * np.pi)
k_x = np.fft.fftshift(k_x)
E_k_x_avg = np.fft.fftshift(E_k_x_avg)

# Plot energy spectrum
plt.figure()
plt.loglog(k_x[Nx//2+1:], E_k_x_avg[Nx//2+1:], label="Computed Spectrum")

# Fit k^(-5/3) scaling in inertial range
inertial_range = (k_x > k_inertial_min) & (k_x < k_inertial_max)
fit_k = k_x[inertial_range]
fit_E = E_k_x_avg[inertial_range]
fit_intercept = np.log(fit_E[-1]) + 5/3 * np.log(fit_k[-1])
k_ref = np.linspace(fit_k[0], fit_k[-1], 100)
fit_line = np.exp(fit_intercept) * k_ref**(-5/3)

plt.loglog(k_ref, fit_line, '--', color='red', label='$k^{-5/3}$ Fit')

plt.xlabel("Wavenumber k")
plt.ylabel("Energy Spectrum E(k)")
plt.legend()
plt.title("1D Energy Spectrum with $k^{-5/3}$ Scaling Exponent")
plt.show()

#####################################
# 3. Compute 2D Energy Spectrum E(k) #
#####################################

# Compute 2D Fourier Transforms
u_hat_2d = np.fft.fft2(u)
v_hat_2d = np.fft.fft2(v)

# Compute energy spectrum in 2D
E_k_2d = (np.abs(u_hat_2d) ** 2 + np.abs(v_hat_2d) ** 2) / 2
E_k_2d = np.fft.fftshift(E_k_2d)

# Compute scalar wavenumber k
kx = np.fft.fftfreq(Nx, d=Lx / Nx) * (2 * np.pi)
ky = kx.copy()
kx, ky = np.meshgrid(np.fft.fftshift(kx), np.fft.fftshift(ky))
k_mag = np.sqrt(kx**2 + ky**2)

# Shell-averaging
k_bins = np.arange(1, np.max(k_mag), 1)
E_k_shell_avg = np.array([np.mean(E_k_2d[(k_mag >= k-0.5) & (k_mag < k+0.5)]) for k in k_bins])

# Plot 2D energy spectrum
plt.figure()
plt.loglog(k_bins, E_k_shell_avg, label="2D Energy Spectrum")

# Fit k^(-5/3) scaling in inertial range for 2D
inertial_range_2d = (k_bins > 28) & (k_bins < 41)
fit_k_2d = k_bins[inertial_range_2d]
fit_E_2d = E_k_shell_avg[inertial_range_2d]
fit_slope_2d, fit_intercept_2d = np.polyfit(np.log(fit_k_2d), np.log(fit_E_2d), 1)
plt.loglog(fit_k_2d, np.exp(fit_intercept_2d) * fit_k_2d**fit_slope_2d, '--', label=f"Fit: k^({fit_slope_2d:.2f})", color='red')

# Compare with 1D spectrum
plt.loglog(k_x[Nx//2+1:], E_k_x_avg[Nx//2+1:], label="1D Spectrum", alpha=0.7)

plt.xlabel("Wavenumber k")
plt.ylabel("Energy Spectrum E(k)")
plt.legend()
plt.title("2D Energy Spectrum with Shell Averaging and Line Fit")
plt.show()


########################################
# 4. Compute Velocity Correlation Functions #
########################################

import numpy as np
import matplotlib.pyplot as plt

# Load dataset
file_path = "isotropic1024_slice.npz"
data = np.load(file_path)

# Extract velocity fields
u = data['u']  # Velocity in x-direction
v = data['v']  # Velocity in y-direction

Nx = u.shape[1]
Lx = 2 * np.pi
dx = Lx / Nx

# Compute mean velocity fields
u_mean = np.mean(u)
v_mean = np.mean(v)

# Compute velocity fluctuations
u_prime = u - u_mean
v_prime = v - v_mean

# Initialize correlation functions
longitudinal_corr = np.zeros(Nx // 2)
transverse_corr = np.zeros(Nx // 2)

# Calculate longitudinal and transverse correlation functions
for r in range(1, Nx // 2):
    u_shift = np.roll(u_prime, -r, axis=1)
    v_shift = np.roll(v_prime, -r, axis=1)
    
    # Longitudinal correlation (u component)
    longitudinal_corr[r] = np.mean(u_prime * u_shift) / np.mean(u_prime**2)
    
    # Transverse correlation (v component)
    transverse_corr[r] = np.mean(v_prime * v_shift) / np.mean(v_prime**2)

# Generate separation distances
r_values = np.arange(1, Nx // 2) * dx

# Plot Longitudinal and Transverse Correlation Functions
plt.figure(figsize=(10, 6))
plt.plot(r_values, longitudinal_corr[1:], label='Longitudinal Correlation', color='blue')
plt.plot(r_values, transverse_corr[1:], label='Transverse Correlation', color='red')

# Plot Formatting
plt.xlabel("Separation Distance r")
plt.ylabel("Correlation Function")
plt.legend()
plt.title("Longitudinal and Transverse Velocity Correlation Functions")
plt.grid(True)

# Show plot
plt.show()


###################################
# 5a. Compute Structure Functions S_p(r) #
###################################

def structure_function(u, p_vals, max_r):
    Sp_r = {p: np.zeros(max_r) for p in p_vals}
    for r in range(1, max_r):
        delta_u = np.abs(u[:, r:] - u[:, :-r])
        for p in p_vals:
            Sp_r[p][r] = np.mean(delta_u ** p)
    return Sp_r

p_values = np.arange(1, 8)
Sp_r_values = structure_function(u, p_values, Nx//2)

plt.figure()
for p in p_values:
    plt.loglog(range(1, Nx//2), Sp_r_values[p][1:], label=f"S_{p}(r)")
plt.legend()
plt.title("Structure Functions S_p(r)")
plt.show()

#################################
# 5b. Verify Kolmogorov's 4/5 Law #
#################################

S3_computed = Sp_r_values[2][1:]

# Step 1: Compute u_rms
u_rms = np.sqrt(np.mean(u**2))

# Given length scale value
length_scale = 1.364

# Step 2: Calculate epsilon using u_rms^3 / length_scale
epsilon = u_rms**3 / length_scale

# print(f"u_rms: {u_rms}")
# print(f"Length Scale: {length_scale}")
# print(f"Epsilon: {epsilon}")

Lx = 2 * np.pi
Nx = u.shape[1]
r_values = np.arange(1, Nx // 2)

# Step 1: Compute u_rms
u_rms = np.sqrt(np.mean(u**2))

# Given length scale value
length_scale = 1.364

# Step 2: Calculate epsilon using u_rms^3 / length_scale
epsilon = u_rms**3 / length_scale

S3_theoretical = -(4 / 5) * epsilon * r_values * (Lx / Nx)

# Compute the third-order structure function S3(r) (Assuming pre-computed S3_r)
def structure_function(u, p, max_r):
    Sp_r = np.zeros(max_r)
    for r in range(1, max_r):
        delta_u = np.abs(u[:, r:] - u[:, :-r])
        Sp_r[r] = np.mean(delta_u ** p)
    return Sp_r

S3_r = structure_function(u, 3, Nx//2)

#############################################
# Step 4: Plot the results
plt.figure()
plt.loglog(r_values * (Lx / Nx), np.abs(S3_r[1:]), label="Computed $S_3(r)$", color='blue')
plt.loglog(r_values * (Lx / Nx), np.abs(S3_theoretical), '--', label="4/5 Law $S_3(r)$", color='red')

# Highlighting the inertial range
k_min, k_max = 28, 47
plt.loglog(r_values[k_min:k_max] * (Lx / Nx), np.abs(S3_r[k_min:k_max]), color='orange', linewidth=3, label='Inertial Range')

# Labeling
plt.xlabel("Separation Distance r")
plt.ylabel("$S_3(r)$")
plt.legend()
plt.title("Verification of Kolmogorov's 4/5 Law")
plt.grid(True)

# Show plot
plt.show()


#################################
# 5c. Extended Self-Similarity (ESS) #
#################################

plt.figure()
for p in p_values:
    if p != 3:
        plt.loglog(S3_computed, Sp_r_values[p][1:], label=f"S_{p}(r) vs S_3(r)")
plt.legend()
plt.title("Extended Self-Similarity (ESS)")
plt.show()

#################################
# 5d. Scaling exponents and error bars #
#################################

# Structure function computation using longitudinal structure
def structure_function(u, p_vals, max_r):
    Sp_r = {p: np.zeros(max_r) for p in p_vals}
    for r in range(1, max_r):
        delta_u = np.abs(u[:, r:] - u[:, :-r])
        for p in p_vals:
            Sp_r[p][r] = np.mean(delta_u ** p)
    return Sp_r

# Define p values and maximum r
p_values = np.arange(1, 8)
max_r = Nx // 2

# Calculate structure functions
Sp_r_values = structure_function(u, p_values, max_r)
S3_r = Sp_r_values[3][1:]

# Perform line fitting for the scaling range in ESS
scaling_range = (S3_r > 1e-3) & (S3_r < 1e-1)
fit_results = {}

plt.figure()
for p in p_values:
    Sp_r_p = Sp_r_values[p][1:]
    plt.loglog(S3_r, Sp_r_p, label=f"$S_{p}(r)$ vs $S_3(r)$")
    
    # Perform linear fit on log-log data
    fit_log_S3 = np.log(S3_r[scaling_range])
    fit_log_Sp = np.log(Sp_r_p[scaling_range])
    slope, intercept = np.polyfit(fit_log_S3, fit_log_Sp, 1)
    fit_results[p] = (slope, intercept)
    plt.loglog(S3_r[scaling_range], np.exp(intercept) * S3_r[scaling_range]**slope, '--')

plt.xlabel("$S_3(r)$")
plt.ylabel("$S_p(r)$")
plt.legend()
plt.title("Extended Self-Similarity (ESS) with Line Fits")
plt.grid(True)
plt.show()

# Extract scaling exponents and calculate errors
zeta_p = []
error_bars = []
for p in p_values:
    slope, _ = fit_results[p]
    zeta_p.append(slope)
    error_bars.append(np.abs(slope - (p / 3)))


# Plot scaling exponents with error bars and connect the points
plt.figure()
plt.errorbar(p_values, zeta_p, yerr=error_bars, fmt='o', label='Computed $zeta_p$', color='yellow')
plt.plot(p_values, zeta_p, color='yellow', linestyle='-', marker='o')  # Connect points
plt.plot(p_values, p_values / 3, 'r--', label="Kolmogorov's Prediction ($zeta_p = p/3$)")

# Labeling
plt.xlabel('$p$')
plt.ylabel('$zeta_p$')
plt.legend()
plt.title("Scaling Exponents $zeta_p$ vs p with Error Bars")
plt.grid(True)

# Display the plot
plt.show()

