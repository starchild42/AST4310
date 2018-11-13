 # importing useful libraries (you may need more)
import numpy as np # numerical package
import matplotlib.pyplot as plt # plotting package
from matplotlib import rc
rc('font',**{'family':'serif'}) # This is for Latex writing

# definition of constants and variables
k = 8.61734e-5 # Boltzmann constant
# add all constants here so you do not need to put them in every function


'''
# definition of functions
def name_function(parameters):
    statements
    statements
    return output_parameters

def another_function(another_params):
    statements
    return another_output

# Main part calling the functions to do things
name_function(parameters)
another_function(another_params)
'''


import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline



chiion = np.array([7, 16, 31, 51]) # Schadee ionization energies into numpy array
k = 8.61734e-5 # Boltzmann constant in eV/deg
temp = 5000. # the decimal point here also makes it a float
u = np.zeros(4) # declare a 4 zero-element array

for r in range(4):
    for s in range(chiion[r]):
        u[r] = u[r] + np.exp(-s / k / temp)
        print(u) # prints all the values of u array (now is not zeros)




def partfunc_E(temp):
    chiion = np.array([7, 16, 31, 51]) # Schadee ionization energies into numpy array
    k = 8.61734e-5 # Boltzmann constant in eV/deg
    u = np.zeros(4) # declare a 4 zero-element array
    for r in range(4):
        for s in range(chiion[r]):
            u[r] = u[r] + np.exp(-s / k / temp)
    return u # returns all the values of u array

# Notice that the variable temp is not inside the function since it will be called when calling
# the function using the command




def boltz_E(temp, r, s):
    u = partfunc_E(temp)
    KeV = 8.61734e-5 # This constant does need to be defined here again if it was before
    relnrs = 1. / u[r - 1] * np.exp(-(s - 1) / (KeV * temp))
    return relnrs



for s in range(1,11): # now the loop starts at 1 and finishes at 10
    print (boltz_E(5000., 1., s))



def saha_E(temp, elpress, ionstage):
    kerg = 1.380658e-16
    kev =8.61734e-5
    h = 6.62607e-27 #erg/s
    elmass = 9.109390e-28 #gram
    kevT = kev * temp
    kergT = kerg * temp
    eldens = elpress / kergT
    chiion = np.array([7, 16, 31, 51 ])

    u = partfunc_E(temp)
    u = np.append(u, 2) # With this command we are adding a new element to the array
    sahaconst = (2. * np.pi * elmass * kergT / (h**2))**1.5 * 2. / eldens
    nstage = np.zeros(5)
    nstage[0] = 1. # We set the first element of the array to a value 1
    
    for r in range(4):
        nstage[r + 1] = nstage[r] * sahaconst * u[r + 1] / u[r] * np.exp(-chiion[r] / kevT)
        ntotal = np.sum(nstage)
        nstagerel = nstage / ntotal
    return nstagerel[ionstage - 1]



temp= 5000.


partfunc_E(temp)
# So, the variable temp has to be defined before executing this command (outside the function)


plt.plot(temp)


for r in range(1,6):
    print (saha_E(20000,1e3,r))


for r in range(1,6):
    print ((saha_E(20000,1e1,r)))



