import numpy as np
from matplotlib import pyplot as plt

# Load in data file
data = np.loadtxt("datafiles/bindingenergies.dat")

# Make arrays containing x-axis and binding energies as function of 
x = data[:,2]
bexpt = data[:,3]
liquiddrop = data[:,4]

# Make plots
plt.plot(x, bexpt ,'b-o', x, liquiddrop, 'r-o')
plt.axis([0,270,-1, 10.0])
plt.xlabel(r'$A$')
plt.ylabel(r'Binding energies in [MeV]') 
plt.legend(('Experiment','Liquid Drop'), loc='upper right') 
plt.title(r'Binding energies from experiment and liquid drop') 
plt.savefig('bindingenergies.pdf') 
plt.savefig('bindingenergies.png')
plt.show()