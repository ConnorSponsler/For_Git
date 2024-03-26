import numpy as np

def calc_differential_flux(cumulative_flux):
  """
  Calculates the differential flux (rate of change of cumulative flux)
  from a provided cumulative flux data structure.

  Args:
      cumulative_flux (numpy.ndarray): A 4-dimensional array containing
          cumulative neutrino flux data. The dimensions should be ordered as
          coszen, height, energy, and flavor (in that order).

  Returns:
      numpy.ndarray: A 4-dimensional array containing the differential flux
          calculated using the average flux within each layer approach.
  """

  # Extract dimensions
  num_coszen, num_heights, num_energies, num_flavors = cumulative_flux.shape

  # Initialize differential flux array
  differential_flux = np.zeros_like(cumulative_flux)

  # Loop through coszen angles
  for i_coszen in range(num_coszen):
    # Loop through energy bins
    for i_energy in range(num_energies - 1):
      # Loop through neutrino flavors
      for i_flavor in range(num_flavors):
        # Calculate layer thickness (assuming constant thickness)
        layer_thickness = np.diff(cumulative_flux[i_coszen, :, i_energy, i_flavor])[0]

        # Calculate average cumulative flux within each layer
        avg_cumulative_flux = (cumulative_flux[i_coszen, 1:, i_energy, i_flavor] +
                              cumulative_flux[i_coszen, :-1, i_energy, i_flavor]) / 2

        # Calculate differential flux (assuming constant within the layer)
        differential_flux[i_coszen, 1:-1, i_energy, i_flavor] = (
            avg_cumulative_flux - cumulative_flux[i_coszen, :-1, i_energy, i_flavor]) / layer_thickness

  return differential_flux



def load_atm_flux_data(filename):
  """
  Loads atmospheric particle flux data from a text file.

  Args:
      filename (str): The path to the text file containing the data.

  Returns:
      numpy.ndarray: A 4-dimensional array containing the cumulative flux data.
          The dimensions should be ordered as coszen, height, energy, and flavor 
          (in that order).
  """

  # Read data from the text file
  data = np.loadtxt(filename, delimiter=' ')  # Assuming comma-separated values

  # Reshape data into a 4D array based on your data structure
  # Modify the reshape parameters based on your actual data organization
  num_coszen, num_heights, num_energy_bins, num_flavors = data.shape[:4]
  return data.reshape(num_coszen, num_heights, num_energy_bins, num_flavors)

# Load your data from the text file
cumulative_flux = load_atm_flux_data("ATM_Particle_Flux_January.txt")

# Calculate differential flux
differential_flux = calc_differential_flux(cumulative_flux)

# Write differential flux to a file
np.savetxt("diff_flux.txt", differential_flux)

