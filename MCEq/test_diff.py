from MCEq.core import config, MCEqRun
import sys
import numpy as np
import crflux.models as crf
# matplotlib used plotting. Not required to run the code.
import matplotlib.pyplot as plt

default_stdout = sys.stdout

months = ['January']#, 'February','March','April','May','June','July','August','September','October','November','December']
h_pts = 100
coszen_pts = 51 
coszen_grid = np.linspace(-1, 1, coszen_pts)
zen_grid = np.arccos(coszen_grid)
zen_deg_grid = np.degrees(zen_grid)

h_grid = np.logspace(6.778, 2,h_pts) # altitudes from 60km to 1m  (in cm)
for month in months:
    f = open("Full_Flux_"+month+".txt", "w")
    mceq = MCEqRun( interaction_model='SIBYLL23C', primary_model = (crf.HillasGaisser2012, 'H3a'), theta_deg = zen_deg_grid[0], density_model=('MSIS00_IC',('SouthPole',month)),
    )
    
    for izen in range(coszen_pts):
        mceq.set_theta_deg(zen_deg_grid[izen])
        
        X_grid = mceq.density_model.h2X(h_grid)
        
        mceq.solve(int_grid=X_grid)
        
            # Loop over the energy grid
        for ien, energy in enumerate(mceq.e_grid):
            # Calculate the flux for all altitudes at this energy
            nue_flux = []
            numu_flux = []
            nutau_flux = []
            antinue_flux = []
            antinumu_flux = []
            antinutau_flux = []
            for idx in range(h_pts):
                nue_flux.append(mceq.get_solution('nue', grid_idx=idx, mag=0)[ien])
                numu_flux.append(mceq.get_solution('numu', grid_idx=idx, mag=0)[ien])
                nutau_flux.append(mceq.get_solution('nutau', grid_idx=idx, mag=0)[ien])
                antinue_flux.append(mceq.get_solution('antinue', grid_idx=idx, mag=0)[ien])
                antinumu_flux.append(mceq.get_solution('antinumu', grid_idx=idx, mag=0)[ien])
                antinutau_flux.append(mceq.get_solution('antinutau', grid_idx=idx, mag=0)[ien])
            nue_flux = np.array(nue_flux)
            numu_flux = np.array(numu_flux)
            nutau_flux = np.array(nutau_flux)
            antinue_flux = np.array(antinue_flux)
            antinumu_flux = np.array(antinumu_flux)
            antinutau_flux = np.array(antinutau_flux)
        
            # Write the energy and differential flux at each altitude to the file
            for idx, altitude in enumerate(h_grid):
                f.write(f"{coszen_grid[izen]} {h_grid[idx]} {mceq.e_grid[ien]} {nue_flux[idx]} {antinue_flux[idx]} {numu_flux[idx]} {antinumu_flux[idx]} {nutau_flux[idx]} {antinutau_flux[idx]}\n")
        
f.close()
