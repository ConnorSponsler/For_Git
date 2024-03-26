import matplotlib.pyplot as plt
import numpy as np

def load_data(filename):
    """
    Loads data from a file. Returns energy, coszen, and separate arrays for
    each flavor's flux.
    """
    data = np.loadtxt(filename)
    if data.shape[1] == 5:  # First file format
        coszen = data[:, 0]
        energy = data[:, 1]
        fluxes = data[:, 2:5]  # Extracting nu_e, nu_mu, nu_tau fluxes
    else:  # Second file format (assuming it has 4 columns)
        # For the second file, we're interested in the last three columns only
        fluxes = data[:, 1:4]  # Adjusted to extract columns 1, 2, 3
        coszen = None
        energy = None
    return coszen, energy, fluxes

# Load data from both files
filename1 = 'fluxes_from_60.txt'  # First data file
filename2 = 'fluxes_from_20.txt'  # Second data file

# Load and filter the data from the first file
coszen1, energy1, fluxes1 = load_data(filename1)
filter_mask = coszen1 == 0  # Apply filter, e.g., for coszen = 0
filtered_energy = energy1[filter_mask]
filtered_fluxes1 = fluxes1[filter_mask]

# Load data from the second file
_, _, fluxes2 = load_data(filename2)

# Ensure the flux arrays from the first file and the entire second file have compatible shapes for division
assert filtered_fluxes1.shape[0] == fluxes2.shape[0], "Filtered data from file 1 and data from file 2 have mismatched number of rows."

# Compute the ratios of the filtered fluxes from the first file to the corresponding fluxes from the second file
flux_ratios = filtered_fluxes1-fluxes2  # Element-wise division

# Plotting the ratios
plt.figure(figsize=(8, 6))
plt.semilogx(filtered_energy, flux_ratios[:, 0], label="Nu_e ratio")
plt.semilogx(filtered_energy, flux_ratios[:, 1], label="Nu_mu ratio")
plt.semilogx(filtered_energy, flux_ratios[:, 2], label="Nu_tau ratio")
plt.xlabel("Energy (GeV)", labelpad=10)
plt.ylabel("Flux Ratio", labelpad=10)
#plt.yscale('log')
plt.title("Neutrino Flux Ratio vs Energy (Log Scale) - Cos(zen) = 0")
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()

# Show the plot
plt.show()
