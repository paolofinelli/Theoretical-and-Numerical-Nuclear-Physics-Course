import numpy as np
from matplotlib import pyplot as plt

# Load in data file
data = np.loadtxt("datafiles/snox.dat")

# Make arrays containing x-axis and binding energies as function of 
x = data[:,1]
y = data[:,3]

# Make plots
plt.plot(x, y,'b-+',markersize=6)
plt.axis([4,18,-7, 12.0])
plt.xlabel(r'Number of neutrons $N$',fontsize=20)
plt.ylabel(r'$\Delta S_n$ [MeV]',fontsize=20)
plt.legend(('Shell gap energies for oxygen isotpes'), loc='upper right')
plt.title(r'Shell gap energies for the oxygen isotopes')
plt.savefig('gapoxygen.pdf')
plt.savefig('gapoxygen.png')
plt.show()
