import matplotlib.pyplot as plt
import numpy as np

def load_data(filename):
  """
  Loads data from a file with energy in the first column and other columns
  representing different flavors. Returns energy, and separate arrays for
  each flavor's flux.
  """
  data = np.loadtxt(filename)
  energy = data[:, 0]
  flux_nu_e = data[:, 1]
  flux_nu_mu = data[:, 2]
  flux_nu_tau = data[:, 3]
  return energy, flux_nu_e, flux_nu_mu, flux_nu_tau

# Replace 'your_data_file.txt' with the actual filename
filename = 'fluxes_from_20.txt'

# Load data from the file
energy, flux_nu_e, flux_nu_mu, flux_nu_tau = load_data(filename)

# Total flux (sum of all flavors)
total_flux = flux_nu_e + flux_nu_mu + flux_nu_tau

# Normalize each flavor flux by the total flux
flux_nu_e_normalized = flux_nu_e / total_flux
flux_nu_mu_normalized = flux_nu_mu / total_flux
flux_nu_tau_normalized = flux_nu_tau / total_flux

# Configure plot
plt.figure(figsize=(8, 6))
plt.semilogx(energy, flux_nu_e_normalized, label="Nu_e")
plt.semilogx(energy, flux_nu_mu_normalized, label="Nu_mu")
plt.semilogx(energy, flux_nu_tau_normalized, label="Nu_tau")
plt.xlabel("Energy (GeV)", labelpad=10)
plt.ylabel("Normalized Flux", labelpad=10)
plt.yscale('log')
plt.title("Normalized Neutrino Flux vs Energy (Log Scale)")
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()

# Show the plot
plt.show()

