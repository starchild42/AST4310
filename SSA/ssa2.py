from __future__ import division #prevents integer division
from pylab import *             #everything we need
from scipy import special

#constants, cgs units
k = 8.61734e-5      #boltzmann
m_e = 9.10939e-28   #electron mass
m_p = 1.67262e-24   #proton mass
m_A = 1.66054e-24   #atomic mass unit
a_0 = 0.529178e-9   #first bohr radius
xi_H_eV = 13.598    #hydrogen ionization energy
xi_H_erg = xi_H_eV*1.60219e-12 #converting to cgs


k_erg = 1.380658e-16    #Boltzmann constant, ergs
h = 6.62607e-27     #Planck constant, ergs
c = 2.99792e10      #light velocity, cgs
temp = 5000. # the decimal point here also makes it a float
KeV = 8.61734e-5


#2.4:

def partitionE(temp, k):
	chiion = array([7, 16, 31, 51]) # Schadee ionization energies into numpy array
	u = zeros(4) # declare a 4 zero-element array
	for r in range(4):
		for s in range(chiion[r]):
			u[r] = u[r] + exp(-s / k / temp)
	return u


def boltz_E(temp, r, s):
	u = partitionE(temp,k)
	relnrs = 1. / u[r - 1] * exp(-(s - 1) / (KeV * temp))
	return relnrs


'''
for s in range(1,11):
	print(boltz_E(5000., 1., s))

0.90181507784
0.0885447147611
0.00869376295073
0.000853597128269
8.38104353107e-05
8.22892771584e-06
8.07957280039e-07
7.93292867443e-08
7.78894613717e-09
7.64757688081e-10
'''



def saha_E(temp, elpress, ionstage):
	kerg = 1.380658e-16	#Boltzmann constant, ergs
	KeVT = KeV * temp
	kergT = kerg * temp
	eldens = elpress / kergT
	
	chiion = array([7, 16, 31, 51 ])
	#chiion = array([6.113, 11.871, 50.91, 67.15])

	u = partitionE(temp,k)
	u = append(u, 2)

	sahaconst = (2. * pi * m_e * kergT / (h**2))**1.5 * 2. / eldens
	nstage = zeros(5)
	nstage[0] = 1.

	for r in range(4):
		nstage[r + 1] = nstage[r] * sahaconst * u[r + 1] / u[r] * exp(-chiion[r] / KeVT)
	ntotal = sum(nstage)
	nstagerel = nstage / ntotal

	return nstagerel[ionstage - 1]



for r in range(1,6):
	print (saha_E(20000,1e3,r))
for r in range(1,6):
	print(saha_E(20000,1e1,r))


'''
2.72775113242e-10
0.000180278462885
0.632005363273
0.36781263824
1.71975186581e-06

7.28751533469e-16
4.81635603475e-08
0.0168847836656
0.982655715769
0.000459452401834
'''


#2.5:
def sahabolt_E(temp, elpress, ion, level):
	return saha_E(temp, elpress, ion) * boltz_E(temp, ion, level)


'''
for s in range(1,6):
	print (sahabolt_E(5000,1e3,1,s))
for s in range(1,6):
	print (sahabolt_E(20000,1e3,1,s))
for s in range(1,6):
	print (sahabolt_E(10000,1e3,2,s))
for s in range(1,6):
	print (sahabolt_E(20000,1e3,4,s))



0.817093536375
0.0802263300862
0.00787702233903
0.000773405450093
7.59368152695e-05
1.22187492385e-10
6.83971553597e-11
3.82868227344e-11
2.14318970926e-11
1.19969791218e-11
0.648954237872
0.203346475021
0.0637175728123
0.0199655739529
0.00625610998149
0.161921366802
0.0906390716843
0.0507372280636
0.0284012872566
0.0158982496407
'''




temp = np.arange(0,30001,1000)
population = np.zeros((5,31))
for T in arange(1,31):
	for r in arange(1,5):
		population[r,T] = sahabolt_E(temp[T],131,r,1) #Change s for energy level
colors = ['-b', '-y', '-g', '-r']
labellst = ['ground stage', 'first ion stage', 'second ion stage', 'third ion stage']
figure(0)



# ground state plot nr 1
figure(2)
for i in range(1, 5):
    plot(temp, population[i, :], label=labellst[i - 1])

xlabel('temperature', size=14)
ylabel('population', size=14)
yscale('log')
ylim([1e-3, 1.1])
legend(loc='best')
title('Population vs Temperature s=1')
show()



#ground-state plot
for i in range(1,5):
	plot(temp,population[i,:],colors[i-1], label=labellst[i-1])


for T in arange(1,31):
	for r in arange(1,5):
		population[r,T] = sahabolt_E(temp[T],131,r,2) #Change s for energy level

#First excited plot
for i in range(1,5):
	plot(temp,population[i,:],colors[i-1])


for T in arange(1,31):
	for r in arange(1,5):
		population[r,T] = sahabolt_E(temp[T],131,r,4) #Change s for energy level
		

#4th plot
for i in range(1,5):
	plot(temp,population[i,:],colors[i-1])
	

yscale('log')
title('Population plot')
ylim([1e-3,1.1])
legend(loc='lower right')
xlabel('temperature',size=14)
ylabel('population',size=14)
savefig('sahabolt_E_variation')
show()






#2.7:

def sahabolt_H(temp,elpress,level):
	keVT = KeV*temp
	kergT = k_erg*temp
	eldens = elpress/kergT

	# energy levels and weights for hydrogen
	nrlevels = 100 # reasonable partition function cut-off value
	g = zeros((2,nrlevels)) # declarations weights (too many for proton)
	chiexc = zeros((2,nrlevels)) # declaration excitation energies (idem)

	for s in range(nrlevels):
		g[0,s] = 2.*(s+1.)**2. # statistical weights
		chiexc[0,s] = 13.598*(1.-1./(s+1.)**2.) # excitation weights
	g[1,0] = 1. # statistical weights free proton
	chiexc[1,0] = 0.

	# partition functions
	u = zeros([2])
	
	for s in range(nrlevels):
		u[0] = u[0] + g[0,s]*exp(-chiexc[0,s]/keVT)
	u[1] = g[1,0]

	# Saha
	sahaconst = (2*np.pi*m_e*kergT /(h*h))**(1.5)*2./eldens
	nstage = zeros(2)
	nstage[0] = 1.
	nstage[1] = nstage[0] * sahaconst * u[1]/u[0] * exp(-13.598/keVT)
	ntotal = sum(nstage) # sum both stages = total hydrogen density

	# Boltzmann
	nlevel = nstage[0]*g[0,level-1]/u[0]*exp(-chiexc[0,level-1]/keVT)
	nlevelrel = nlevel/ntotal # fraction of total hydrogen density
	
	return nlevelrel


print (sahabolt_H(6000,1e2,1))
'''
for s in range(6):
	print( s+1, g[0,s], chiexc[0,s], g[0,s]*exp(-chiexc[0,s]/keVT))

for s in range(0,nrlevels,10):
	print( s+1, g[0,s], chiexc[0,s], g[0,s]*exp(-chiexc[0,s]/keVT))

'''



#2.8 Solar Ca+K versus Ha: line strength

temp = arange(1000,20001,100)
CaH = zeros(temp.shape)
Caabund = 2.0e-6

for i in range(0,191):
	NCa = sahabolt_E(temp[i],1e2,2,1) # is equal to sahabolt_Ca
	NH = sahabolt_H(temp[i],1e2,2)
	CaH[i] = NCa*Caabund/NH

plot(temp,CaH, label=r'strength ratio Ca$^+$K / H$\alpha$')
yscale('log')
title('Strength Ratio')
xlabel(r'temperature $T / K$', size=14)
ylabel(r'Ca II K / H$\alpha$', size=14)
legend(fontsize=14)
show()

print ('Ca/H ratio at 5000 K = ', CaH[np.argwhere(temp==5000)][0][0])


#2.9 9 Solar Ca+K versus Ha: temperature sensitivity

temp = arange(2000,12001,100)
dNCadT = zeros(temp.shape)
dNHdT = zeros(temp.shape)
dT = 1.

for i in range(101):
	NCa = sahabolt_E(temp[i],1e2,2,1)
	NCa2 = sahabolt_E(temp[i]-dT,1e2,2,1)
	dNCadT[i] = (NCa - NCa2)/(dT*NCa)
	NH = sahabolt_H(temp[i],1e2,2)
	NH2 = sahabolt_H(temp[i]-dT,1e2,2)
	dNHdT[i] = (NH-NH2)/(dT*NH)


figure()
plot(temp,absolute(dNHdT), label=r'H')
plot(temp,absolute(dNCadT), label=r'Ca$^+$K')
title('Temperature Sensitivity')
yscale('log')
#ylim(1e-9,1)
xlabel(r'temperature $T/K$', size=14)
ylabel(r'$\left| \left( \Delta n(r,s) / \Delta T \right) / n(r,s) \right|$', size=20)
legend(loc=4, fontsize=12)


NCa = zeros(temp.shape)
NH = zeros(temp.shape)

for i in range(101):
	NCa[i] = sahabolt_E(temp[i],1e2,2,1)
	NH[i] = sahabolt_H(temp[i],1e2,2)


plot(temp,NH/np.amax(NH), ls='--', label = 'rel. population. H')
plot(temp,NCa/np.amax(NCa), ls='--', label = r'rel. population. Ca$^+$')
title('Relative Populations')
show()



#2.10 Hot stars versus cool stars

for T in arange(2e3,2e4+1,2e3):
	print (T, sahabolt_H(T,1e2,1))

temp = arange(1e3,2e4+1,1e2)
nH = zeros(temp.shape)

for i in range(191):
	nH[i] = sahabolt_H(temp[i],1e2,1)


plot(temp,nH)
title('H Fraction for Temperature')
xlabel('temperature $T/K$', size=14)
ylabel('neutral hydrogen fraction', size=14)

legend()
show()


