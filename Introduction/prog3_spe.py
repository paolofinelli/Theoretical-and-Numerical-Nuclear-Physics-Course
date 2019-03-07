import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

# Load in data file
data = np.loadtxt("datafiles/bindingenergies.dat")

# Make arrays containing A on x-axis and binding energies 
A = data[:,2]
bexpt = data[:,3]

# The function we want to fit to, only two terms here
def func(A,a1, a2): return a1*A-a2*(A**(2.0/3.0))

# function to perform nonlinear least square with guess for a1 and a2 
popt, pcov = curve_fit(func, A, bexpt, p0 = (16.0, 18.0))
a1 = popt[0]
a2 = popt[1]
liquiddrop = a1*A-a2*(A**(2.0/3.0))

plt.plot(A, bexpt ,'bo', A, liquiddrop, 'ro')
plt.axis([0,270,-1, 10.0])
plt.xlabel(r'$A$')
plt.ylabel(r'Binding energies in [MeV]')
plt.legend(('Experiment','Liquid Drop'), loc='upper right')
plt.title(r'Binding energies from experiment and liquid drop')
plt.savefig('bindingenergies.pdf')
plt.savefig('bindingenergies.png')
plt.show()