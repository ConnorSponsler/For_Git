from MCEq.core import config, MCEqRun
import sys
import numpy as np
import crflux.models as crf
# matplotlib used plotting. Not required to run the code.
import matplotlib.pyplot as plt

default_stdout = sys.stdout

months = ['January']
h_pts = 100
coszen_pts = 50 
coszen_grid = np.linspace(-1, 1, coszen_pts)
zen_grid = np.arccos(coszen_grid)
zen_deg_grid = np.degrees(zen_grid)

h_grid = np.logspace(6.778, 2,h_pts) # altitudes from 60km to 1m  (in cm)
for month in months:
    f = open("ATM_Particle_Flux_"+month+".txt", "w")
    mceq = MCEqRun( interaction_model='SIBYLL23C', primary_model = (crf.HillasGaisser2012, 'H3a'), theta_deg = zen_deg_grid[0], density_model=('MSIS00_IC',('SouthPole',month)),
    )
    
    for izen in range(coszen_pts):
        mceq.set_theta_deg(zen_deg_grid[izen])
        
        X_grid = mceq.density_model.h2X(h_grid)
        
        mceq.solve(int_grid=X_grid)
        
        for idx in range(h_pts):
            numu_flux = (mceq.get_solution('numu', grid_idx=idx))
            nue_flux = (mceq.get_solution('nue', grid_idx=idx))
            for ien in range(len(mceq.e_grid)):
                sys.stdout = f
                print(coszen_grid[izen], ' ', h_grid[idx], ' ', mceq.e_grid[ien], ' ', nue_flux[ien], ' ', numu_flux[ien] )
                sys.stdout = default_stdout
        
