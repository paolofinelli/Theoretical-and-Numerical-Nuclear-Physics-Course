import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.units as units
import matplotlib.ticker as ticker
#rc('text',usetex=True) 
#rc('font',**{'family':'serif','serif':['Woods-Saxon potential']}) 
#font = {'family' : 'serif',
#        'color'  : 'darkred',
#        'weight' : 'normal',
#        'size'   : 16,
#        }
        
v0 = 50
A = 100
a = 0.5
r0 = 1.25
R = r0*A**(0.3333)
x = np.linspace(0.0, 10.0)
y = -v0/(1+np.exp((x-R)/a))

plt.plot(x, y, 'b-')
plt.title(r'Woods-Saxon potential', fontsize=20)
plt.text(3, -40, r'Parameters: $A=20$, $V_0=50$ [MeV]')
plt.text(3, -44, r'$a=0.5$ [fm], $r_0=1.25$ [fm]')
plt.xlabel(r'$r$ [fm]',fontsize=20)
plt.ylabel(r'$V(r)$ [MeV]',fontsize=20)

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.savefig('woodsaxon.pdf', format='pdf')
plt.show()