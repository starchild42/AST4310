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




#3.1
def planck(temp, wav):
    B_wav = (2*h*c**2/(wav)**5)*( 1/(exp( (h*c)/(wav*k_erg*temp) ) - 1) )
    return B_wav



wav = arange(1000, 20801, 200)
b = zeros(wav.shape)

xlabel(r'wavelength $\lambda / \AA$', size=14)
ylabel(r'Planck function', size=14)
xlim(0, 20800)

for T in range(8000, 5000-1, -200):
    b[:] = planck(T, wav[:]*1e-8)
    plot(wav,b, label='%g Kelvin' %T)
legend(loc='best')
savefig('planck.png')
show()

#3.2:
B = 2
tau = arange(0.01, 10.01, 0.01)
intensity = zeros(tau.shape)

for I0 in range(4, -1, -1):
    intensity[:] = I0*exp(-tau[:]) + B*(1-exp(-tau[:]))
    plot(tau, intensity, label = 'intensity I0 = ' + str(I0))

xlabel(r'optical depth $\tau$', size=14)
ylabel('intensity', size=14)
legend(fontsize=12)
savefig('depth.png')
show()


#3.3

def voigt(gamma, x):
    z = (x+1j*gamma)
    V = special.wofz(z).real
    return V

#voigt profile:
u = np.arange(-10, 10.1, 0.1)
a = array([0.001, 0.01, 0.1, 1])
vau = zeros((a.shape[0], u.shape[0]))

for i in range(4):
    vau[i,:] = voigt(a[i],u[:])
    plot(u[:], vau[i,:], label = 'a = ' + str(a[i]))

ylim(0,1)
xlim(-10,10)
legend(fontsize = 12)
ylabel('voigt profile', size=12)
savefig('voigt.png')
show()

for i in range(4):
    vau[i,:] = voigt(a[i],u[:])
    plot(u[:], vau[i,:], label = 'a = ' + str(a[i]))

yscale('log')
legend(fontsize=12, loc = 8)
xlabel('u', size=14)
ylabel('logarithmic voigt profile', size=12)
savefig('voigt_log.png')
show()


#schuster-schwarzchild line profile

Ts = 5700
T1 = 4200
a = 0.1
wav = 5000e-8
tau0 = 1
u = arange(-10,10.1,0.1)
intensity = zeros(u.shape)

for i in range(201):
    tau = tau0 * voigt(a, u[i])
    intensity[i] = planck(Ts, wav) * exp(-tau) + planck(T1, wav) * (1-exp(-tau))

plot(u, intensity)
xlabel('u')
ylabel('Schuster-Schwarzschild profile')
savefig('sch_sch.png')
show()

logtau0 = arange(-2, 2.1, 0.5)

for itau in range(9):
    for i in range(201):
        tau = 10**(logtau0[itau]) * voigt(a, u[i])
        intensity[i] = planck(Ts, wav) * exp(-tau) + planck(T1,wav)*(1-exp(-tau))
    plot(u, intensity, label = r'$\log{(\tau_0)} = $' + np.str(logtau0[itau]))

legend(loc=3, fontsize=12)
xlabel('u')
ylabel('Schuster-Schwarzschild profile')
savefig('sch_sch_var.png')
show()

for iwav in range(1,4):
    wav = (iwav**2+1)*1e-5
    for itau in range(8):
        for i in range(201):
            tau = 10**(logtau0[itau]) * voigt(a, u[i])
            intensity[i] = planck(Ts, wav) * exp(-tau) + planck(T1, wav)*(1-exp(-tau))
        intensity = intensity/intensity[0]
        plot(u, intensity[:], linewidth=1)

xlabel('u')
ylabel('Schuster-Schwarzschild profile')
savefig('sch_sch_obs.png')
show()



#3.4
def profile(a, tau0, u):
    Ts = 5700
    T1 = 4200
    wav = 5000e-8
    usize = u.size
    intensity = zeros(usize)
    for i in range(usize):
        tau = tau0 * voigt(a, u[i])
        intensity[i] = planck(Ts, wav)*exp(-tau) + planck(T1, wav)*(1-exp(-tau))
    return intensity




u = arange(-200, 200.4, 0.4)
a = 0.1
tau0 = 1e2
intensity = profile(a, tau0, u)

plot(u, intensity)
xlabel('dimensionless wavelength u')
ylabel('intensity')
savefig('intensity_e2.png')
show()

reldepth = (intensity[0] - intensity)/(intensity[0])
plot(u, reldepth)
xlabel('dimensionless wavelength u')
ylabel('intensity')
savefig('reldepth_e2.png')
show()
eqw = sum(reldepth)*0.4
print (eqw)






#3.5:
u = arange(-200, 200.4, 0.4)
a = 0.1

tau0 = logspace(-2, 4, 61)
eqw = zeros(tau0.size)

for i in range(61):
    intensity = profile(a, tau0[i], u)
    reldepth = (intensity[0] - intensity)/intensity[0]
    eqw[i] = sum(reldepth)*0.4
    print (eqw[i])

plot(tau0, eqw)
xlabel(r'$\tau_0$', size=18)
ylabel(r'equivalent width $W_{\lambda}$', size=14)
xscale('log')
yscale('log')
savefig('width.png')
show()
