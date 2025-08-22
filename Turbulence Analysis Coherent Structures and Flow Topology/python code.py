# ========================================
# Coherent Structures and Flow Topology
# ========================================

#ans 2  PDF of λ1, λ2 and λ3 -----------------------------------------
# Plot PDFs of eigenvalues λ1, λ2, λ3

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

# Load velocity data from file
data = np.load('isotropic1024_stack3.npz')
ux = data['u']  
uy = data['v']
uz = data['w']

# Extract middle plane (z=1) for each velocity component
ux_mid = ux[:, :, 1]
uy_mid = uy[:, :, 1]
uz_mid = uz[:, :, 1]

# Grid parameters
L = 2 * np.pi
Nx = 1024
dx = dy = dz = L / Nx

# Compute x-derivatives (periodic boundaries)
ux_x = (np.roll(ux_mid, -1, axis=0) - np.roll(ux_mid, 1, axis=0)) / (2 * dx)
uy_x = (np.roll(uy_mid, -1, axis=0) - np.roll(uy_mid, 1, axis=0)) / (2 * dx)
uz_x = (np.roll(uz_mid, -1, axis=0) - np.roll(uz_mid, 1, axis=0)) / (2 * dx)

# Compute y-derivatives (periodic boundaries)
ux_y = (np.roll(ux_mid, -1, axis=1) - np.roll(ux_mid, 1, axis=1)) / (2 * dy)
uy_y = (np.roll(uy_mid, -1, axis=1) - np.roll(uy_mid, 1, axis=1)) / (2 * dy)
uz_y = (np.roll(uz_mid, -1, axis=1) - np.roll(uz_mid, 1, axis=1)) / (2 * dy)

# Compute z-derivatives using central differencing (non-periodic)
ux_z = (ux[:, :, 2] - ux[:, :, 0]) / (2 * dz)
uy_z = (uy[:, :, 2] - uy[:, :, 0]) / (2 * dz)
uz_z = (uz[:, :, 2] - uz[:, :, 0]) / (2 * dz)

# Construct velocity gradient tensor A_ij for each point (3x3)
A = np.empty((Nx, Nx, 3, 3))
A[..., 0, 0] = ux_x  # dux/dx
A[..., 0, 1] = ux_y  # dux/dy
A[..., 0, 2] = ux_z  # dux/dz
A[..., 1, 0] = uy_x  # duy/dx
A[..., 1, 1] = uy_y  # duy/dy
A[..., 1, 2] = uy_z  # duy/dz
A[..., 2, 0] = uz_x  # duz/dx
A[..., 2, 1] = uz_y  # duz/dy
A[..., 2, 2] = uz_z  # duz/dz

# Compute eigenvalues for each gradient tensor
eigenvalues = np.linalg.eigvals(A)

# Handling complex eigenvalues 
if np.iscomplexobj(eigenvalues):
    #eigenvalues = np.real(eigenvalues)
    eigenvalues = np.abs(eigenvalues)
    #print("Note: Some eigenvalues were complex; real parts taken.")
    print("Note: Some eigenvalues were complex; absolute value is taken.")

# Sort eigenvalues in descending order (λ1 > λ2 > λ3)
eigenvalues_sorted = np.sort(eigenvalues, axis=2)[:, :, ::-1]


# Flatten eigenvalues for PDF calculation
lambda1 = eigenvalues_sorted[:, :, 0].flatten()
lambda2 = eigenvalues_sorted[:, :, 1].flatten()
lambda3 = eigenvalues_sorted[:, :, 2].flatten()

# Create common bins for all eigenvalues
all_lambdas = np.concatenate([lambda1, lambda2, lambda3])
global_min, global_max = np.nanmin(all_lambdas), np.nanmax(all_lambdas)
bins = np.linspace(global_min, global_max, 100)
bin_centers = (bins[:-1] + bins[1:]) / 2

# Compute histograms (density normalized)
hist1, _ = np.histogram(lambda1, bins=bins, density=True)
hist2, _ = np.histogram(lambda2, bins=bins, density=True)
hist3, _ = np.histogram(lambda3, bins=bins, density=True)

# Avoid log(0) by replacing zeros with NaN
hist1[hist1 == 0] = np.nan
hist2[hist2 == 0] = np.nan
hist3[hist3 == 0] = np.nan

# Plot PDFs on semilog scale
plt.figure(figsize=(10, 6))
plt.semilogy(bin_centers, hist1, 'r-', label=r'$\lambda_1$')
plt.semilogy(bin_centers, hist2, 'g-', label=r'$\lambda_2$')
plt.semilogy(bin_centers, hist3, 'b-', label=r'$\lambda_3$')

plt.xlabel('Eigenvalue', fontsize=12)
plt.ylabel('Probability Density (log scale)', fontsize=12)
plt.title('PDFs of Eigenvalues of Velocity Gradient Tensor', fontsize=14)
plt.legend()
plt.grid(True, which='both', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()

#########################################################################
#ans 3  -----------------------------------------
# Verify whether ⟨|P|⟩ ≈ 0 over the 2D domain

# Compute invariants P, Q, R
#P = np.trace(A, axis1=2, axis2=3)  # P = λ1 + λ2 + λ3
P= -(ux_x + uy_y + uz_z)  # P = λ1 + λ2 + λ3


# Question 3: Verify incompressibility (P ≈ 0)
mean_abs_P = np.mean((np.abs(P)))
print(f"Mean |P|: {mean_abs_P}")

mean_P = np.mean((P))
print(f"Mean P: {mean_P}")


#########################################################################
#ans 4 -----------------------------------------

# Plot Q/Q_rms and R/R_rms

Q = -0.5 * np.trace(A @ A, axis1=2, axis2=3)  # Q = -0.5 * trace(A^2)
R = -1/3 * np.trace(A @ A @ A, axis1=2, axis2=3)  # R = -1/3 * trace(A^3), if a+b+c=0 then a^3+b^3+c^3=3abc


# #comparing Q with different methods
# # Compute S_ij and Omega_ij
# S = 0.5 * (A + np.transpose(A, axes=(0, 1, 3, 2)))  # Symmetric part
# Omega = 0.5 * (A - np.transpose(A, axes=(0, 1, 3, 2)))  # Antisymmetric part

# # Compute Q from S and Omega
# Q_alt = 0.5 * (np.sum(Omega**2, axis=(2, 3)) - np.sum(S**2, axis=(2, 3)))

# # Check consistency with earlier Q
# print("Difference in Q calculations:", np.max(np.abs(Q - Q_alt)))
# print('Q_alt',np.mean(Q_alt))
# print('Q',np.mean(Q))

Q_rms = np.sqrt(np.mean(Q**2))
R_rms = np.sqrt(np.mean(R**2))

plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.imshow(Q / Q_rms, cmap='RdBu_r', vmin=-5, vmax=5)
plt.colorbar(label=r'$Q / Q_{\mathrm{rms}}$')
plt.title('Normalized Q Field')

plt.subplot(1, 2, 2)
plt.imshow(R / R_rms, cmap='RdBu_r', vmin=-5, vmax=5)
plt.colorbar(label=r'$R / R_{\mathrm{rms}}$')
plt.title('Normalized R Field')
plt.tight_layout()
plt.show()

#########################################################################
# #ans 5 -----------------------------------------

# Joint distribution of Q and R
# Compute enstrophy ω = ∇ × u
omega_x = uz_y - uy_z
omega_y = ux_z - uz_x
omega_z = uy_x - ux_y
enstrophy = omega_x**2 + omega_y**2 + omega_z**2
Q_w_mean = np.mean(enstrophy) / 4  # <Q_w> = <ω²>/4

# Normalize Q and R
Q_norm = Q / Q_w_mean
R_norm = R / (Q_w_mean**1.5)

# Create 2D histogram (logarithmic bins)
H, xedges, yedges = np.histogram2d(
    Q_norm.flatten(), R_norm.flatten(),
    bins=100, range=[[-10, 10], [-10, 10]], density=True
)
H = gaussian_filter(H, sigma=1)  # Smooth for visualization
H[H == 0] = np.nan  # Mask zeros for log scale

# Plot joint distribution
plt.figure(figsize=(8, 6))
plt.contourf(xedges[:-1], yedges[:-1], np.log10(H.T), levels=20, cmap='viridis')
plt.colorbar(label='log10(Probability Density)')

# Plot discriminant line 27R² + 4Q³ = 0
x = np.linspace(-10, 0, 100)
plt.plot(x, -2 * np.sqrt(-x**3 / 27), 'r--', label=r'$27R^2 + 4Q^3 = 0$')
plt.plot(x, 2 * np.sqrt(-x**3 / 27), 'r--')

plt.xlabel(r'$Q / \langle Q_w \rangle$')
plt.ylabel(r'$R / \langle Q_w \rangle^{3/2}$')
plt.title('Joint Distribution in Q-R Plane')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)
plt.xlim(-10, 10)
plt.ylim(-10, 10)
plt.tight_layout()
plt.show()


#########################################################################

