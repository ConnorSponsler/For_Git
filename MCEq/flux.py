from MCEq.core import MCEqRun
import numpy as np
import crflux.models as crf

# Define the altitude grid (in cm)
h_pts = 100
h_grid = np.linspace(6000000, 100, h_pts)  # Altitudes from 60km to 1m

# Initialize MCEqRun object with specified parameters
mceq = MCEqRun(interaction_model='SIBYLL23C',
               primary_model=(crf.HillasGaisser2012, 'H3a'),
               theta_deg=0,
               density_model=('MSIS00_IC', ('SouthPole', 'January')))

# Convert altitude grid to depth grid
X_grid = mceq.density_model.h2X(h_grid)

# Solve the cascade equations for the entire grid
mceq.solve(int_grid=X_grid)

# Open a file to save the differential flux data
with open("nuProd_January.txt", "w") as f:
    # Loop over the energy grid
    for ien, energy in enumerate(mceq.e_grid):
        # Calculate the flux for all altitudes at this energy
        nue_flux_at_energy = []
        numu_flux_at_energy = []
        nutau_flux_at_energy = []
        for idx in range(h_pts):
            nue_flux_at_energy.append(mceq.get_solution('nue', grid_idx=idx, mag=0)[ien])
            numu_flux_at_energy.append(mceq.get_solution('numu', grid_idx=idx, mag=0)[ien])
            nutau_flux_at_energy.append(mceq.get_solution('nutau', grid_idx=idx, mag=0)[ien])
        nue_flux_at_energy = np.array(nue_flux_at_energy)
        numu_flux_at_energy = np.array(numu_flux_at_energy)
        nutau_flux_at_energy = np.array(nutau_flux_at_energy)

        # Calculate the differential flux with respect to height
        #(negative because cumulative production (flux) increases as height decreases)
        diff_nue_flux_wrt_height = -1*np.gradient(nue_flux_at_energy, h_grid)
        diff_numu_flux_wrt_height = -1*np.gradient(numu_flux_at_energy, h_grid)
        diff_nutau_flux_wrt_height = -1*np.gradient(nutau_flux_at_energy, h_grid)
        
        # Write the energy and differential flux at each altitude to the file
        for idx, altitude in enumerate(h_grid):
            f.write(f"{altitude/100000} {energy} {diff_nue_flux_wrt_height[idx]} {diff_numu_flux_wrt_height[idx]} {diff_nutau_flux_wrt_height[idx]}\n")
            
            
# Open a file to save the differential flux data
with open("nuFlux_January.txt", "w") as f:
    # Loop over the energy grid
    for ien, energy in enumerate(mceq.e_grid):
        # Calculate the flux for all altitudes at this energy
        nue_flux_at_energy = []
        numu_flux_at_energy = []
        nutau_flux_at_energy = []
        for idx in range(h_pts):
            nue_flux_at_energy.append(mceq.get_solution('nue', grid_idx=idx, mag=0)[ien])
            numu_flux_at_energy.append(mceq.get_solution('numu', grid_idx=idx, mag=0)[ien])
            nutau_flux_at_energy.append(mceq.get_solution('nutau', grid_idx=idx, mag=0)[ien])
        nue_flux_at_energy = np.array(nue_flux_at_energy)
        numu_flux_at_energy = np.array(numu_flux_at_energy)
        nutau_flux_at_energy = np.array(nutau_flux_at_energy)
        
        # Write the energy and differential flux at each altitude to the file
        for idx, altitude in enumerate(h_grid):
            f.write(f"{altitude/100000} {energy} {nue_flux_at_energy[idx]} {numu_flux_at_energy[idx]} {nutau_flux_at_energy[idx]}\n")
            