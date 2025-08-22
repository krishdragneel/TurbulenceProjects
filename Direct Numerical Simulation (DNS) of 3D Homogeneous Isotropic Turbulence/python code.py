import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

def compute_turbulence_parameters(u, v, w, nu_kinematic, L, dx):
    """Compute Reynolds number and Kolmogorov scales."""
    U_rms = np.sqrt(np.mean(u**2))
    Re = (U_rms * L) / nu_kinematic
    eta = (nu_kinematic**3 / (U_rms**3 / L))**0.25  # Kolmogorov length scale
    nu_eta = (nu_kinematic * (U_rms**3 / L))**0.5  # Kolmogorov velocity scale
    tau_eta = (nu_kinematic / (U_rms**3 / L))**0.5  # Kolmogorov time scale
    grid_ratio = dx / eta  # Grid cell size compared to Kolmogorov length scale
    
    print(f"Reynolds Number (Re): {Re:.2f}")
    print(f"Kolmogorov Length Scale (eta): {eta:.6f}")
    print(f"Kolmogorov Velocity Scale (u_eta): {nu_eta:.6f}")
    print(f"Kolmogorov Time Scale (tau_eta): {tau_eta:.6f}")
    print(f"Grid Cell Size / Kolmogorov Length: {grid_ratio:.2f}")

def plot_normalized_kinetic_energy_field(u, v, w, Nx, Ny):
    """Plot the 2D normalized kinetic energy field."""
    k = 0.5 * (u**2 + v**2 )
    k_avg = np.mean(k)  # Spatial average of k
    w_field = k / k_avg  # Normalized kinetic energy field
    
    plt.figure(figsize=(8, 6))
    plt.imshow(w_field, origin='lower', cmap='jet', extent=[0, Nx, 0, Ny])
    plt.colorbar(label='Normalized Kinetic Energy')
    plt.title('2D Normalized Kinetic Energy Field')
    plt.xlabel('Grid Number (x)')
    plt.ylabel('Grid Number (y)')
    plt.show()

def plot_normalized_vorticity_field(u, v, dx, Nx, Ny):
    """Compute and plot the normalized vorticity field."""
    dvdx = (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1)) / (2 * dx)
    dudy = (np.roll(u, -1, axis=0) - np.roll(u, 1, axis=0)) / (2 * dx)
    omega_z = dvdx - dudy  # Vorticity component
    omega_rms = np.sqrt(np.mean(omega_z**2))  # Root-mean-square vorticity
    omega_z_normalized = omega_z / omega_rms  # Normalized vorticity
    
    plt.figure(figsize=(8, 6))
    plt.imshow(omega_z_normalized, origin='lower', cmap='seismic', extent=[0, Nx, 0, Ny])
    plt.colorbar(label='Normalized Vorticity')
    plt.title('2D Normalized Vorticity Field')
    plt.xlabel('Grid Number (x)')
    plt.ylabel('Grid Number (y)')
    plt.show()

def plot_pdfs(u, v, w, omega_z, dx):
    """Plot PDFs for given quantities."""
    quantities = {
        #'PDF of Velocity Components (u, v, w)': [u, v, w],
        
        'PDF of Vorticity Component (ω_z)': [omega_z],
        'PDF of Velocity Gradients (∂u/∂x)': [(np.roll(u, -1, axis=1) - np.roll(u, 1, axis=1)) / (2 * dx),
                                                    ],
        'PDF of Velocity Gradients (∂v/∂y)': [
                                                    (np.roll(v, -1, axis=0) - np.roll(v, 1, axis=0)) / (2 * dx)],                                            
        'PDF of Enstrophy (ω_z^2)': [omega_z**2]
    }
    
    for title, data_list in quantities.items():
        plt.figure(figsize=(8, 6))
        for data in data_list:
            data_std = (data - np.mean(data)) / np.std(data)
            pdf, bins = np.histogram(data_std.flatten(), bins=100, density=True)
            plt.plot(bins[:-1], pdf, label=title)
        
        x = np.linspace(-4, 4, 100)
        plt.plot(x, stats.norm.pdf(x, 0, 1), 'k--', label='N(0,1)')
        plt.yscale('log')
        plt.xlabel('Normalized Value')
        plt.ylabel('Probability Density')
        plt.legend()
        plt.title(title)
        plt.show()

def plot_pdfs_velocities(u, v, w, omega_z, dx):
    """Plot PDFs for velocity components u, v, w in a single graph."""
    plt.figure(figsize=(8, 6))
    labels = ['u', 'v', 'w']
    colors = ['r', 'g', 'b']
    
    for data, label, color in zip([u, v, w], labels, colors):
        data_std = (data - np.mean(data)) / np.std(data)
        pdf, bins = np.histogram(data_std.flatten(), bins=100, density=True)
        plt.plot(bins[:-1], pdf, label=label, color=color)
    
    x = np.linspace(-4, 4, 100)
    plt.plot(x, stats.norm.pdf(x, 0, 1), 'k--', label='N(0,1)')
    plt.yscale('log')
    plt.xlabel('Normalized Value')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.title('PDF of Velocity Components (u, v, w)')
    plt.show()



def plot_cdfs(u, v, w, omega_z, dx):
    """Plot CDFs for given quantities."""
    quantities = {
        'CDF of Velocity Components (u, v, w)': [u, v, w],
        'CDF of Vorticity Component (ω_z)': [omega_z],
        'CDF of Velocity Gradients (∂u/∂x, ∂v/∂y)': [
            (np.roll(u, -1, axis=1) - np.roll(u, 1, axis=1)) / (2 * dx),
            (np.roll(v, -1, axis=0) - np.roll(v, 1, axis=0)) / (2 * dx)
        ],
        'CDF of Enstrophy (ω_z^2)': [omega_z**2]
    }
    
    for title, data_list in quantities.items():
        plt.figure(figsize=(8, 6))
        if title == 'CDF of Velocity Components (u, v, w)':
            labels = ['u', 'v', 'w']
            colors = ['r', 'g', 'b']
            for data, label, color in zip(data_list, labels, colors):
                data_std = (data - np.mean(data)) / np.std(data)
                sorted_data = np.sort(data_std.flatten())
                cdf = np.arange(len(sorted_data)) / float(len(sorted_data))
                plt.plot(sorted_data, cdf, label=label, color=color)
        
        elif title == 'CDF of Velocity Gradients (∂u/∂x, ∂v/∂y)':
            labels = ['∂u/∂x', '∂v/∂y']
            colors = ['m', 'c']
            for data, label, color in zip(data_list, labels, colors):
                data_std = (data - np.mean(data)) / np.std(data)
                sorted_data = np.sort(data_std.flatten())
                cdf = np.arange(len(sorted_data)) / float(len(sorted_data))
                plt.plot(sorted_data, cdf, label=label, color=color)
        
        else:
            for data in data_list:
                data_std = (data - np.mean(data)) / np.std(data)
                sorted_data = np.sort(data_std.flatten())
                cdf = np.arange(len(sorted_data)) / float(len(sorted_data))
                plt.plot(sorted_data, cdf, label=title)
        
        plt.xlabel('Normalized Value')
        plt.ylabel('Cumulative Probability')
        plt.legend()
        plt.title(title)
        plt.show()

def compute_and_print_skewness_kurtosis(u, v, w, omega_z):
    """Compute and print skewness and kurtosis for given quantities."""
    dudx = (np.roll(u, -1, axis=1) - np.roll(u, 1, axis=1)) / (2 * dx)
    dvdy = (np.roll(v, -1, axis=0) - np.roll(v, 1, axis=0)) / (2 * dx)
    quantities = {
        #'Velocity Components (u, v, w)': [u, v, w],
        'Velocity Components (u)': [u],
        'Velocity Components (v)': [v],
        'Velocity Components (w)': [w],
        'Vorticity Component (ω_z)': [omega_z],
        'Enstrophy Component (ω_z*w_z)': [omega_z**2],
        'Velocity Gradient (∂u/∂x)': [dudx],
        'Velocity Gradient (∂v/∂y)': [dvdy],
    }
    
    for label, data_list in quantities.items():
        for data in data_list:
            data_std = (data - np.mean(data)) / np.std(data)
            skewness = stats.skew(data_std.flatten())
            kurtosis = stats.kurtosis(data_std.flatten())
            print(f"{label} Skewness: {skewness:.6f}")
            print(f"{label} Kurtosis: {kurtosis:.6f}")

# Load data
data = np.load('isotropic1024_slice.npz')

# Extract velocity components
u = data['u']  # x-component of velocity
v = data['v']  # y-component of velocity
w = data['w']  # z-component of velocity







def plot_pdf_by_x_ranges(u, x_ranges, Ny):
    """Plot PDFs of u(x, y) for different x ranges."""
    plt.figure(figsize=(8, 6))
    for x_start, x_end in x_ranges:
        u_subset = u[x_start:x_end, :].flatten()
        u_std = (u_subset - np.mean(u_subset)) / np.std(u_subset)
        pdf, bins = np.histogram(u_std, bins=100, density=True)
        plt.plot(bins[:-1], pdf, label=f'x ∈ [{x_start}, {x_end}]')
    
    x = np.linspace(-4, 4, 100)
    plt.plot(x, stats.norm.pdf(x, 0, 1), 'k--', label='N(0,1)')
    plt.yscale('log')
    plt.xlabel('Normalized Value')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.title('PDF of u(x, y) for Different x Ranges')
    plt.show()

def plot_local_reynolds_number_field(u, v, nu, dx, Nx, Ny):
    """Compute and plot the local Reynolds number field."""
    dudx = (np.roll(u, -1, axis=1) - np.roll(u, 1, axis=1)) / (2 * dx)
    dudy = (np.roll(u, -1, axis=0) - np.roll(u, 1, axis=0)) / (2 * dx)
    d2udx2 = (np.roll(u, -1, axis=1) - 2*u + np.roll(u, 1, axis=1)) / (dx**2)
    d2udy2 = (np.roll(u, -1, axis=0) - 2*u + np.roll(u, 1, axis=0)) / (dx**2)
    
    Nlin = u * dudx + v * dudy  # Nonlinear term
    Visc = nu * (d2udx2 + d2udy2)  # Viscous term
    
    Re_local = np.abs(Nlin) / np.abs(Visc + 1e-10)  # Local Reynolds number , and avoid dividing by 0
    
    

    skip = 5  # Adjust this to control memory usage
    fields = [
        
        (Re_local, 'Local Reynolds Number Field', 'jet', 0, 10)
    ]

    for field, title, cmap, vmin, vmax in fields:
        plt.figure(figsize=(8, 6))
        plt.pcolor(field[::skip, ::skip], cmap=cmap, vmin=vmin, vmax=vmax)
        plt.colorbar(label=title)
        plt.title(title)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig(f"{title.replace(' ', '_')}.png", bbox_inches='tight', dpi=240)
        plt.show()




# Given parameters
nu_kinematic = 0.000185  # Kinematic viscosity
L = 1.364  # Integral length scale
Nx = Ny = 1024  # Grid size
Lx = Ly = 2 * np.pi  # Domain size
dx = Lx / Nx  # Grid cell size

# Compute vorticity
omega_z = (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1)) / (2 * dx) - (np.roll(u, -1, axis=0) - np.roll(u, 1, axis=0)) / (2 * dx)

# Compute and print turbulence parameters
compute_turbulence_parameters(u, v, w, nu_kinematic, L, dx)

# # Plot normalized kinetic energy field
plot_normalized_kinetic_energy_field(u, v, w, Nx, Ny)

# # Plot the normalized vorticity field
plot_normalized_vorticity_field(u, v, dx, Nx, Ny)

# # Plot PDFs
plot_pdfs_velocities(u, v, w, omega_z, dx)
plot_pdfs(u, v, w, omega_z, dx)

# # Plot CDFs
plot_cdfs(u, v, w, omega_z, dx)

# # Plot PDFs of u(x, y) for different x ranges
plot_pdf_by_x_ranges(u, [(0, 128), (129, 256), (257, 512)], Ny)

# # Plot Local Reynolds Number field
plot_local_reynolds_number_field(u, v, nu_kinematic, dx, Nx, Ny)

# Compute and print skewness and kurtosis
compute_and_print_skewness_kurtosis(u, v, w, omega_z)

