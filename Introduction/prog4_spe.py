import numpy as np
from matplotlib import pyplot as plt

# Load in data file
data = np.loadtxt("datafiles/snox.dat")

# Make arrays containing x-axis and binding energies as function of 
x = data[:,1]
y = data[:,2]

# Make plots
plt.plot(x, y,'b-+',markersize=6)
plt.axis([4,18,-1, 25.0])
plt.xlabel(r'Number of neutrons $N$',fontsize=20)
plt.ylabel(r'$S_n$ [MeV]',fontsize=20)
plt.legend(('Separation energies for oxygen isotpes'), loc='upper right')
plt.title(r'Separation energy for the oxygen isotopes')
plt.savefig('snoxygen.pdf')
plt.savefig('snoxygen.png')
plt.show()
