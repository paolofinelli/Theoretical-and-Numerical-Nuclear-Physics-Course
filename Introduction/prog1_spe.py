import numpy as np
from matplotlib import pyplot as plt

# Load in data file
data = np.loadtxt("datafiles/bindingenergies.dat")

# Make arrays containing x-axis and binding energies as function of A 
x = data[:,2]
bexpt = data[:,3]

# Make plots
plt.plot(x, bexpt ,'ro')
plt.axis([0,270,-1, 10.0])
plt.xlabel(r'$A$')
plt.ylabel(r'Binding energies in [MeV]')
plt.legend(('Experiment'), loc='upper right')
plt.title(r'Binding energies from experiment')
plt.savefig('expbindingenergies.pdf')
plt.savefig('expbindingenergies.png')
plt.show()