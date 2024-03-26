from MCEq.core import MCEqRun
import numpy as np
import crflux.models as crf

# Initialize MCEqRun object with specified parameters
mceq = MCEqRun(interaction_model='SIBYLL23C',
               primary_model=(crf.HillasGaisser2012, 'H3a'),
               theta_deg=0,
               density_model=('MSIS00_IC', ('SouthPole', 'January')))

# Solve the cascade equations for the entire grid
mceq.solve()

#save surface flux
with open("th0_cumu-flux.txt", "w") as f:
    # Loop over the energy grid
    for ien, energy in enumerate(mceq.e_grid):

        nue_flux = mceq.get_solution('nue', mag=0)[ien]
        numu_flux = mceq.get_solution('numu', mag=0)[ien]
        nutau_flux = mceq.get_solution('nutau', mag=0)[ien]
        antinue_flux = mceq.get_solution('antinue', mag=0)[ien]
        antinumu_flux = mceq.get_solution('antinumu', mag=0)[ien]
        antinutau_flux = mceq.get_solution('antinutau', mag=0)[ien]
        
        f.write(f"{energy} {nue_flux} {numu_flux} {nutau_flux} {antinue_flux} {antinumu_flux} {antinutau_flux}\n")
            
            
