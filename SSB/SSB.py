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
kerg = 1.380658E-16	# Boltzmann constant in erg K
kJoule = kerg*1e-7	# Boltzmann constant in joule/K
R_s = 696300  		# Radius sun in km
D_s = 1.496e8		# Distance sun-earth in km

#1.1: plot height and temp
plt.plot(h, temp)
plt.xlabel('height [km]')
plt.ylabel('temperature [k]')
plt.axis([-500, 2500,2000,10000])
#plt.show()



# 1.2: plot total pressure and column mass linearly
plt.plot(colm, ptot)
plt.xlabel('Column mass')
plt.ylabel('p total')
#plt.show()


# Logarithmic
plt.plot(colm,ptot)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('log(Column mass)')
plt.ylabel('log(p total)')
#plt.show()


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
#plt.show()

# Remaining metal fractions:
mfrac = 1. - average((rho_h+rho_he)/dens)
print 'Remaining Metal Fraction = ', mfrac



#Plot of Column mass against height:
plt.plot(h, colm)
plt.title('Column Mass against Height')
plt.ylabel(r'm [g cm$^{-2}$]')
plt.xlabel(r'Height [km]')
#plt.show()

#Logarithmic:
plt.plot(h, colm)
plt.yscale('log')
plt.title('Column Mass against Height')
plt.ylabel(r'log(m) [g cm$^{-2}$]')
plt.xlabel(r'Height [km]')
#plt.show()


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
#plt.show()


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



#-----------------EARTH-----------------------------------

# reading earth.dat
Eh, ElogP, Etemp, Elogdens, ElogN = loadtxt("earth.dat",usecols=(0,1,2,3,4), unpack=True)


#plot temperature, pressure, density and gas density

plt.plot(Eh, ElogP)  #ratio curves
plt.title('Air pressure vs height logarithmic')
plt.ylabel(r'log$P$ [dyn cm$^{-2}$]')
plt.xlabel(r'Height [km]')
#plt.show()


plt.plot(Eh, Etemp)  #ratio curves
plt.title('Temperature against height')
plt.ylabel(r'Temperature [K]')
plt.xlabel(r'Height [km]')
#plt.show()


plt.plot(Eh, Elogdens)  #ratio curves
plt.title('Gas density against height logarithmic')
plt.ylabel(r'log $\rho$ [g cm$^{-3}$]')
plt.xlabel(r'Height [km]')
#plt.show()


plt.plot(Eh, ElogN)  #ratio curves
plt.title('Particle density against height')
plt.ylabel(r'log $N$ [cm$^{-3}$]')
plt.xlabel(r'Height [km]')
#plt.show()


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
'''
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

#atmospheric column mass at earth's surface,
#compared to the sun's column mass:
Ecolm =Eptot[where(Eh==0)][0]/g_E
print('col_m earth h=0', Ecolm)
print('col_m sun h=0',colm[where(h==0)][0])

#comparing sun and earth values:
#ratio of particle densities at h=0
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



#energy flux of sunshine reaching our planet:
#Flux= pi*B*T_eff
#iRad = (4*pi*R**2)/(4*pi*D**2)*Flux
#N_phot = (pi*R**2/D**2)*N_phot_top

n_p_s = pi*R_s**2/D_s**2*(20*5770**3/(2*pi))  #photons from the sun
n_p_e = 20*Etemp[where(Eh==0)][0]**3  #photons from Earth
print('sun phot:',n_p_s, 'earth phot:', n_p_e)
print('photon ratio sun/earth:',n_p_s/n_p_e)
print(N[where(Eh==0)][0]) #particle dens surface E
print(n_p_s/N[where(Eh==0)][0]) #ratio between solar photons&particles in E atmosphere

'''



#---------------------------------------------------------

#EXERCISE 2: CONTINUOUS SPECTRUM FROM THE SOLAR ATMOSPHERE

#2.1:
# reading solspect.dat
wav, F,Fcont,I,Icont= loadtxt("solspect.dat",usecols=(0,1,2,3,4), unpack=True)


# to obtain maxima
print('max(Ic)= ',max(Icont),'at',wav[where(Icont == max(Icont))])
# max Ic = 4.59999, at array[0.41]



#Plot four spectral distributions
plt.plot(wav,F)
plt.plot(wav,Fcont)
plt.plot(wav,I)
plt.plot(wav,Icont)
plt.legend([r"$F_{cont}^{\nu}$",r"$F_{cont}^{\nu}$ smoothed",r"$I_{cont}^{\nu}$",r"$F_{cont}^{\nu}$ smoothed"], loc=1)
#plt.xlim(-500,2.5e3)
#plt.ylim(0,10000)
plt.xlabel(r"wavelength $\lambda$ [$\mu m$]")
plt.ylabel(r"intensity, flux [$erg cm^{-2} s^{-1} ster^{-1} Hz^{-1}$]")
plt.title("solar continuum radiation")
plt.show()


#convert to values per frequency
c=3e14
factor = wav**2/c*1e10
print( 'max(Ic) =', max(Icont*factor),'at',wav[where(Icont*factor == max(Icont*factor))])
# max(Ic) = 4.2026666666666675e-05

#Plot four spectral distributions with frequencies
plt.plot(wav,F*factor)
plt.plot(wav,Fcont*factor)
plt.plot(wav,I*factor)
plt.plot(wav,Icont*factor)
plt.legend([r"$F_{cont}^{\nu}$",r"$F_{cont}^{\nu}$ smoothed",r"$I_{cont}^{\nu}$",r"$F_{cont}^{\nu}$ smoothed"], loc=1)
#plt.xlim(-500,2.5e3)
#plt.ylim(0,10000)
plt.xlabel(r"wavelength $\lambda$ [$\mu m$]")
plt.ylabel(r"intensity, flux [$erg cm^{-2} s^{-1} ster^{-1} Hz^{-1}$]")
plt.title("solar continuum radiation, frequency plot")
plt.show()

print('peak:', Icont*factor) #The peak

#Planck function
#Constants in cgs
k = 1.38065e-16
h_P = 6.626076e-27 # Planck constant [erg s]
c = 2.997929e10  # Velocity of light [micro m/s]

def planck(temp,wav):
	blambda = 2*h_P*c**2/(wav**5*(exp(h_P*c/(wav*k*temp))-1))
	return blambda

def brightTemp(wav,I):
	var1= h_P*c/(wav*k)
	var2= 2*h_P*c**2/(I*wav**5)*1e-4
	return var1/log(var2+1)


#plot of brightness temperature:
plt.plot(wav, brightTemp(wav*1e-4,Icont*1e10))
plt.xlabel(r"wavelength $\lambda$ [$\mu m$]")
plt.ylabel(r"$T_b$ [K]")
plt.title("Brightness temperature")
plt.show()



#2.2: Continuous Extinction
h, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot, dens = loadtxt("falc.dat",usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)

'''
#
plt.plot(wav, EXmin)
plt.xlabel(r'Wavelength [$\mu$m]')
plt.ylabel(r'$\kappa /P_e [10^{-26}$cm$^2$/dyn cm$^{-2}$]')
plt.xscale('log')
plt.yscale('log')
plt.show()
'''


def exthmin(wav,temp,eldens):
	# H-minus extinction, from Gray 1992
	# input:
	#  wav = wavelength [Angstrom] (float or float array)
	#  temp = temperature [K]
	#  eldens = electron density [electrons cm-3]
	# output:
	#  H-minus bf+ff extinction [cm^2 per neutral hydrogen atom]
	#  assuming LTE ionization H/H-min
	# physics constants in cgs (all cm)
	k=1.380658e-16	# Boltzmann constant [erg/K]
	h=6.626076e-27	# Planck constant [erg s]
	c=2.997929e10	# velocity of light [cm/s]
	
	# other parameters
	theta=5040./temp
	elpress=eldens*k*temp

	# evaluate H-min bound-free per H-min ion ? Gray (8.11)
	# his alpha = my sigma in NGSB/AFYC (per particle without stimulated)
	sigmabf = (1.99654 -1.18267E-5*wav +2.64243E-6*wav**2-4.40524E-10*wav**3 +3.23992E-14*wav**4-1.39568E-18*wav**5 +2.78701E-23*wav**6)
	sigmabf *= 1e-18 # cm^2 per H-min ion
	if size(wav) > 1:
		sigmabf[where(wav>16444)] = 0 # H-min ionization limit at lambda = 1.6444 micron
	elif (size(wav) == 1):
		if wav> 16444:
			sigmabf = 0

	# convert into bound-free per neutral H atom assuming Saha = Gray p135
	# units: cm2 per neutral H atom in whatever level (whole stage)
	graysaha=4.158E-10*elpress*theta**2.5*10.**(0.754*theta)# Gray (8.12)
	kappabf=sigmabf*graysaha # per neutral H atom
	kappabf=kappabf*(1.-exp(-h*c/(wav*1E-8*k*temp)))# correct stimulated

	# check Gray's Saha-Boltzmann with AFYC (edition 1999) p168
	# logratio=-0.1761-np.log10(elpress)+np.log10(2.)+2.5*np.log10(temp)-theta*0.754
	# print 'Hmin/H ratio=',1/(10.**logratio) # OK, same as Gray factor SB

	# evaluate H-min free-free including stimulated emission = Gray p136
	lwav = log10(wav)
	f0 = - 2.2763 - 1.6850*lwav + 0.76661*lwav**2 - 0.0533464*lwav**3
	f1 =   15.2827 - 9.2846*lwav + 1.99381*lwav**2 - 0.142631*lwav**3
	f2 = - 197.789 + 190.266*lwav - 67.9775*lwav**2 + 10.6913*lwav**3 - 0.625151*lwav**4
	ltheta = log10(theta)
	kappaff = 1e-26*elpress*10**(f0+f1*ltheta+f2*ltheta**2) #Gray(8.13)
	return kappabf+kappaff


#temp = T(h=0), h[69]=0
wav = wav*1e4
nneutH=nhyd-nprot
sigmaT= 6.648e-25
EXmin = exthmin(0.5*1e4,temp,nel)*nneutH 
Thomson = nel*sigmaT
plt.semilogy(h,EXmin)
plt.semilogy(h,Thomson)
plt.semilogy(h,EXmin+Thomson)

plt.legend([r"$H^{-}$",r"Thomson",r"Total"], loc=1)
#plt.xlim(0,2e4)
#plt.ylim(0,2.5e-24)
plt.xlabel(r"Height h [km]")
plt.ylabel(r"\alpha_{\lambda} [$cm^{-1}$]")
plt.title(r"H- extinction at $\lambda = 0.5 \mu m$(FALC.dat)")
plt.show()




#2.3: Optical depth:
tau = zeros(len(tau5), dtype=float) 
nneutH=nhyd-nprot
sigmaT= 6.648e-25
EXmin = exthmin(500,temp,nel) 
Thomson = nel*sigmaT
ext = EXmin+Thomson

for i in range(1,len(tau)):
	tau[i] = tau[i-1] + 0.5*(ext[i]+ext[i-1])*(h[i-1]-h[i])*1e5
# index zero is not accounted for, so tau[0] = 0 because we have already initialized


plt.plot(h,tau5,'--', label = 'tau5')
plt.plot(h,tau, label = 'tau')
plt.yscale('log')
plt.xlabel(r"Height h [km]")
plt.ylabel(r"optical depth[$\tau_{\lambda}$]")
plt.legend([r'$\tau_{5000}$ in FALC',r'$\tau_{5000}$ from H- and Thomson extinction'])
plt.title(r'$\tau_{500}$ from FALC compared to numerical integration')
plt.show()




#2.4: Emergent intensity and height formation:
sigma_Thomson= 6.648e-25 # Thomson cross-section [cm^2]

contfuncCalc=zeros(len(wav))
for j in range(len(wav)):
	ext = zeros(len(tau5))
	tau = zeros(len(tau5))
	integrand = zeros(len(tau5))
	contfunc = zeros(len(tau5))
	intt = 0.0
	hint = 0.0
	for i in range(1, len(tau5)):
		ext[i] = (exthmin(wav[j]*1e4, temp[i], nel[i])*(nhyd[i]-nprot[i])+ sigma_Thomson*nel[i])
		tau[i] = tau[i-1] + 0.5*(ext[i] + ext[i-1])*(h[i-1]-h[i])*1e5
		integrand[i] = planck(temp[i],wav[j]*1e-4)*exp(-tau[i])
		intt += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
		hint += h[i]*0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
		contfunc[i] = integrand[i]*ext[i]
	contfuncCalc[j]=intt
	# note : exthmin has wavelength in [Angstrom], planck in [cm]
	mean = hint / intt
	#tau5[69]=1
plt.plot(wav, contfuncCalc*1e-14)
plt.plot(wav,Icont)

plt.xlabel(r"wavelength $\lambda$  [$\mu$ m]")
plt.ylabel(r"Intensity [$10^{14} erg s^{-1} cm^{-2} ster^{-1} cm^{-1}$]")
plt.title('Observed and computed continuum intensity')
plt.legend([r'computed from FALC',r'observed (Allen 1978)'])
plt.show()
