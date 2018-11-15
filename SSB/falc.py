from numpy import *
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sherif'})


# Reading falc.dat file and making arrays:
h, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot,dens = loadtxt('falc.dat', usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)


# CONSTANTS
m_h = 1.67352e-24
m_he = 3.97*m_h
m_e = 9.10939e-28
m_p = 1.67262e-24
pgas = ptot-dens*vturb**2/2
k = 1.38065e-16
n_phot = 20*temp**3
keV = 8.61734E-5 	# Boltzmann constant in eV/K
kerg = 1.380658E-16	 # Boltzmann constant in erg K
kJoule = kerg*1e-7	# Boltzmann constant in joule/K


#1.1: plot height and temp
plt.plot(h, temp)
plt.xlabel('height [km]')
plt.ylabel('temperature [k]')
plt.axis([-500, 2500,2000,10000])
plt.show()



# 1.2: plot total pressure and column mass linearly
plt.plot(colm, ptot)
plt.xlabel('Column mass')
plt.ylabel('p total')
plt.show()


# Logarithmic
plt.plot(colm,ptot)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('log(Column mass)')
plt.ylabel('log(p total)')
plt.show()


#Find c value on average:
g_S = average(ptot/colm)
print 'g_S = ', g_S


# Plot ratio of hydrogen mass density to total mass density vs height:
rho_h = nhyd*m_h # H mass density
n_he = 0.1*nhyd # He number density
rho_he = n_he*m_he # He mass density

plt.plot(h, rho_h/dens)
plt.plot(h, rho_he/dens)
plt.plot(h, (rho_h+rho_he)/dens)
plt.legend([r'$\rho_{H}/\rho_{total}$',
            r'$\rho_{He}/\rho_{total}$',
            r'$(\rho_{H}+\rho_{He})/\rho_{total}$'], loc='best')
plt.title('Density Fractions VS Height')
plt.ylim([0, 1.1])
plt.ylabel(r'$\rho/\rho_{total}$')
plt.xlabel(r'Height [km]')
plt.show()

# Remaining metal fractions:
mfrac = 1. - average((rho_h+rho_he)/dens)
print 'Remaining Metal Fraction = ', mfrac



#Plot of Column mass against height:
plt.plot(h, colm)
plt.title('Column Mass against Height')
plt.ylabel(r'm [g cm$^{-2}$]')
plt.xlabel(r'Height [km]')
plt.show()

#Logarithmic:
plt.plot(h, colm)
plt.yscale('log')
plt.title('Column Mass against Height')
plt.ylabel(r'log(m) [g cm$^{-2}$]')
plt.xlabel(r'Height [km]')
plt.show()


# Plot gas density against height:
#rho = rho_0*exp(-h/H_rho) #density scale height in the deep photosphere
#first, we need to scale the height:

rhodiv = dens[-1]/exp(1)
rho_e = zeros(size(h))
rho_e.fill(rhodiv)

plt.plot(h,dens)
plt.plot(h,rho_e)
plt.title(r'Gas Density against Height')
plt.xlabel(r'Height [km]')
plt.ylabel(r'Gas density $\rho$ [g cm$^{-3}$]')
plt.legend(['Gas density', r'$\rho/e$'])
plt.show()


#find gas pressure
pgas = pgasptot*ptot #gas pressure
idealgas = (nhyd + nel)*kerg*temp


#plot product
plt.plot(h,pgas*1e-4) #falc
plt.plot(h, idealgas*1e-4) #ideal gas
plt.title(r'Gas Pressure against Height')
plt.legend(['Falc', 'Ideal Gas'])
plt.ylabel(r'Gas pressure [$10^4$ dyn cm$^{-2}$]')
plt.xlabel(r'Height [km]')
#show()


#ratio of curves, plot:
plt.plot(h, pgas/idealgas)
plt.title('Ratio of gas pressures, falc and ideal gas law')
plt.ylabel('Gas pressure ratio')
plt.xlabel('Height [km]')
#show()


#Additional number density:
idealgas_hel = (nhyd + nel + n_he)*kerg*temp

#Do previous plots again with new ideal gas:

#ratio plot:
plt.plot(h, pgas/idealgas_hel)
plt.title('Ratio of gas pressures, falc and ideal gas law')
plt.ylabel('Gas pressure ratio')
plt.xlabel('Height [km]')
plt.ylim([0,1.11])
#show()

#new plot of gas pressures
plt.plot(h,pgas*1e-4) #falc
plt.plot(h, idealgas_hel*1e-4) #ideal gas
plt.title(r'Gas Pressure against Height')
plt.legend(['Falc', 'Ideal Gas w He'])
plt.ylabel(r'Gas pressure [$10^4$ dyn cm$^{-2}$]')
plt.xlabel(r'Height [km]')
#show()


#plot of total H density against height and overplot
#of e density, p density and density of e's not from H ionization:

n_eb = (nhyd-nprot) #bound electron

plt.plot(h,nhyd)
plt.plot(h,nel)
plt.plot(h,nprot)
plt.plot(h,n_eb)
plt.title('Number Densities against Height')
plt.legend([r'$n_H$', r'$n_e$', r'$n_p$', r'$n_{be}$'])
plt.ylabel('n [cm^-3]')
plt.xlabel('Height [km]')
#plt.yscale('log')
#plt.show()


hyd_ion = nprot/nhyd

#plot ionization frac of H, log vs height
plt.plot(h,hyd_ion)
plt.title('Ionization Fraction of Hydrogen Against Height')
plt.ylabel(r'log($n_p/n_H$')
plt.xlabel('Height [km]')
plt.yscale('log')
plt.xlim([-100,2220])
#plt.show()


#photon density calculations
n_phot = 20*temp**3


print ("Photon density at deepest model location: %g" % n_phot[argmin(h)])
print ("Hydrogen density at deepest model location: %g" % nhyd[argmin(h)])
print ("Photon density at highest model location: %g" % (20*5770**3/(2*pi)))
print ("Hydrogen density at highest model location: %g" % nhyd[argmax(h)])



#------------------EARTH------------------------------------------

# reading earth.dat
Eh, ElogP, Etemp, Elogdens, ElogN = loadtxt("earth.dat",usecols=(0,1,2,3,4), unpack=True)


#plot temperature, pressure, density and gas density

plt.plot(Eh, ElogP)  #ratio curves
plt.title('Air pressure vs height logarithmic')
plt.ylabel(r'log$P$ [dyn cm$^{-2}$]')
plt.xlabel(r'Height [km]')
plt.show()


plt.plot(Eh, Etemp)  #ratio curves
plt.title('Temperature against height')
plt.ylabel(r'Temperature [K]')
plt.xlabel(r'Height [km]')
plt.show()


plt.plot(Eh, Elogdens)  #ratio curves
plt.title('Gas density against height logarithmic')
plt.ylabel(r'log $\rho$ [g cm$^{-3}$]')
plt.xlabel(r'Height [km]')
plt.show()


plt.plot(Eh, ElogN)  #ratio curves
plt.title('Particle density against height')
plt.ylabel(r'log $N$ [cm$^{-3}$]')
plt.xlabel(r'Height [km]')
plt.show()


#normalizing units
P = 10.**ElogP
dens = 10.**Elogdens
N = 10.**ElogN
dense = zeros(size(Eh))
dense.fill(dens[where(Eh==0)][0]/exp(1))

#pressure and dens stratification plot all in one
plt.plot(Eh,dens)
plt.plot(Eh,dense)
plt.show()

plt.plot(Eh,dens(max(dens)))
plt.plot(Eh,P(max(P)))
plt.plot(Eh,N/max(N))
plt.yscale('log')
plt.ylabel('Values')
plt.xlabel('Height [km]')
plt.title('Density and pressure stratification')
plt.legend([r'$\rho/\rho_{max}$', '$P/P_{max}$', '$N/N_{max}$'])
plt.show()


#mean molecular weight against height,plot
mu_e=dens/(N*m_h)

plt.plot(Eh,mu_e)
plt.title('Mean Molecular Weight against Height')
plt.xlabel('Height [km]')
plt.ylabel(r'$\mu_e$[]')
plt.show()

#density scale h, lower earth atmosphere
g_E= 980.665 #[cm s**(-2)]
eH_p = kJoule*temp/((mu_e*m_u*1e-3)*(g_E*1e-2))*1e-3
print('scale height lower earth atmosphere:', eH_p[0],'km')


#comparing sun and earth values:
print('------------------------------------------')
print('All Values, h = 0')

print('Matter dens sun  :', dens[where(h == 0)][0])
print('Matter dens earth:', dense[where(Eh == 0)][0])
print('------------------------------------------')

N1 = (nhyd+n_he+n_e)
print('Particle density sun  :', N1[where(h == 0)][0])
print('Particle density earth:', N[where(Eh == 0)][0])
print('------------------------------------------')
Eptot = 10**ElogP

print('Pressure sun  :', ptot[where(h == 0)][0])
print('Pressure earth:', Eptot[where(Eh == 0)][0])
print('------------------------------------------')
print('Temperature sun  :', temp[where(h == 0)][0])
print('Temperature earth:', Etemp[where(Eh == 0)][0])
print('------------------------------------------')
print('Ratio particle densities:', N[where(Eh == 0)][0]/N1[where(h == 0)][0])
print('------------------------------------------')






