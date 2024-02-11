import numpy as np
import matplotlib.pyplot as plt

# Replace 'path_to_your_file.txt' with the actual path to your file
file_path = 'numu_Differential_Flux_January.txt'

# Assuming the file content is space-separated (or adjust the delimiter accordingly)
data = np.loadtxt(file_path, delimiter=' ')

# Separate the columns into individual arrays
heights = data[:, 0]  # Heights
energies = data[:, 1]  # Energies
diff_flux = data[:, 2]  # Differential flux
# Sum of the entire 3rd column (differential flux)
total_flux_sum = np.sum(diff_flux)

# Sums corresponding to each unique height value
unique_heights = np.unique(heights)
height_sums = np.array([np.sum(diff_flux[heights == h]) for h in unique_heights])
densities = 0.0012*np.exp(-(unique_heights)/(7.594*100000))
# Calculate the ratios of the h_i sums over the total
ratios = height_sums / total_flux_sum


fig, ax1 = plt.subplots()

# Plotting the first dataset with height as the y-axis
ax1.plot(ratios, unique_heights, 'o-',color='g')
ax1.set_xlabel('X1 axis',color='g')
ax1.set_ylabel('Height (m)', color='g')

# Create a second x-axis for atmospheric density
ax2 = ax1.twiny()
ax2.plot(densities, unique_heights, 'b-',color='b')
ax2.set_xlabel('Atmospheric Density (kg/m^3)', color='b')

# Format primary x-axis in scientific notation
ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))


plt.show()
