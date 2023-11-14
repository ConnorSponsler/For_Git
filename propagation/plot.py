import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('6_atmos.txt')
names = ['energies', 'Q1','nuQS','Q3', 'nuE', 'nuMu', 'nuTau']
energies = data[:,0]
for i in range(1, data.shape[1]):  # Loop over other columns as y-values
    plt.loglog(energies, data[:, i], label= names[i])


plt.xlabel('Energy (eV)')  # Assuming the energy unit is electron volts (eV)
plt.ylabel('Oscillation Probability')
plt.title('Neutrino Oscillation Probabilities vs Energy')
plt.legend() 
plt.show()
