import numpy as np
from scipy import special



#constants and variables
keV = 8.61734E-5 	# Boltzmann constant in eV/K
kerg = 1.380658E-16	 # Boltzmann constant in erg K
kJoule = kerg*1e-7	 # Boltzmann constant in joule/K
h = 6.62607E-27	 # Planck constant in erg s
c = 2.99792E10  # Speed of light in cm/s
m_e = 9.109390E-28  # Electron mass in grams
m_H = 1.67352E-24   # Hydrogen mass in grams
m_He = 3.98*m_H 	# Helium mass in grams
m_p = 1.67262E-24  # Proton mass in grams
m_u = 1.660538e-24  # Atomic mass unit in grams
m_Na = 22.99*1.6605*1e-24  # Sodium mass in grams
R_sun = 696300  # Radius sun in km
D_sun = 1.496e8 	# Distance sun-earth in km
erg2eV = 6.242e11  # Convert erg to ev




# Use planck func from SSA
# This function takes in wavelength in micrometers and temperature in K, and returns the planck function with units [1e10 erg/cm^2/s/micrometer/steradian]
def planck(temp, wav):
    b = 2.*h*c*c/((wav*1e-4)**5) / (np.exp((h*c) / ((wav*1e-4) * k_erg * temp)) - 1)
    return b*1e-14  # 1e-14 coverts to same units as in tables



# This function takes in wavelength in micrometers and temperature in K, and returns the planck function with units [1e10 erg/cm^2/s/micrometer/steradian]
def planck_nu(temp, wav):
    conversion_factor = 1. / c * (wav * 1e-4) ** 2 * 1e4
    b = 2.*h*c*c/((wav*1e-4)**5) / (np.exp((h*c) / ((wav*1e-4) * k_erg * temp)) - 1)
    return b*1e-14*conversion_factor  # 1e-14 coverts to same units as in tables


# Invert planck function
# Takes in intensity in units 10^10 erg cm^-2 s^-1 ster^-1\mu m^-1, and wavelength in units \mu m, and returns brightness temperature T_b
def planck_invert(intensity, wavelength):
    t_b = h * c / (wavelength * 1e-4 * k_erg) * 1. / \
        np.log((2. * h * c * c) / (intensity * 1e14 * (wavelength * 1e-4) ** 5) + 1.)
    return t_b


# Continuous extinction
def exthmin(wav, temp, eldens):
    # in: 	wav = wavelength [Angstrom] (float or fltarr)
    # temp = temperature [K]
    # eldens = electron density [electrons cm-3]
    #
    # out: 	H-minus bf+ff extinction [cm^2 per neutral hydrogen atom]
    # assuming LTE ionization H/H-min

    theta = 5040. / temp
    elpress = eldens * k_erg * temp
    sigmabf = (1.99654 - 1.18267E-5*wav + 2.64243E-6*wav**2 - 4.40524E-10*wav**3 + 3.23992E-14*wav**4
               - 1.39568E-18*wav**5 + 2.78701E-23*wav**6)

    sigmabf *= 1e-18
    for i in range(np.size(wav)-1):
        if sigmabf[i] >= 16444:
            sigmabf[i] = 0

    graysaha = 4.158E-10*elpress*theta**2.5*10.**(0.754*theta)
    kappabf = sigmabf*graysaha
    kappabf = kappabf*(1. - np.exp(-h * c / (wav * 1E-8 * k_erg * temp)))

    lwav = np.log10(wav)

    f0 = -2.2763 - 1.6850*lwav + 0.76661*lwav**2 - 0.0533464*lwav**3
    f1 = 15.2827 - 9.2846*lwav + 1.99381*lwav**2 - 0.142631*lwav**3
    f2 = (-197.789 + 190.266*lwav - 67.9775*lwav**2 + 10.6913*lwav**3 - 0.625151*lwav**4)

    ltheta = np.log10(theta)
    kappaff = 1E-26*elpress*10**(f0+f1*ltheta+f2*ltheta**2)

    return kappaff + kappabf


# Find nearest function
def find_nearest(array, value):
    index = (np.abs(array-value)).argmin()
    return array[index], index


# Flux integration
def intensity_mu(tau5, wavelength, mu, temp, nel, nhyd, nprot, height):
    # Takes in an array of optical depths tau5, a wavelength and an angle mu and
    # returns the intensity and the height multiplied intensity
    ext = np.zeros(np.size(tau5))
    tau = np.zeros(np.size(tau5))
    integrand = np.zeros(np.size(tau5))
    intt = 0
    hint = 0
    for i in range(1, len(tau)):  # the index zero is not accounted for
        ext[i] = exthmin(wavelength*1e4, temp[i], nel[i])*(nhyd[i]-nprot[i])+0.664e-24*nel[i]
        tau[i] = tau[i-1] + 0.5*(ext[i]+ext[i-1])*(height[i-1]-height[i])*1e5
        integrand[i] = planck(temp[i], wavelength)*np.exp(-tau[i]/mu)/mu
        intt += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
        hint += height[i]*0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
    return [intt, hint]


# Vacuum to air conversion
def vac_to_air(wavelength):
    # Takes in vacuum wavelength in angstrom and returns air wavelength in angstrom
    return 0.99972683*wavelength + 0.0107 - 196.25/wavelength


# partition function sodium
def partfunc_Na(temperature):
    u = np.zeros(3)
    theta = 5040. / temperature
    c0 = 0.30955
    c1 = -0.17778
    c2 = 1.10594
    c3 = -2.42847
    c4 = 1.70721
    logu1 = (c0 + c1 * np.log10(theta) + c2 * np.log10(theta)**2 + c3 * np.log10(theta)**3 + c4 * np.log10(theta)**4)
    u[0] = 10**logu1
    u[1] = 1.
    u[2] = 6.
    return u


# Boltzmann for sodium
def boltzmann_na(temp, r, s):
    # Boltzmann distribution n_r,s/N_r
    e_n1 = h*c / 5895.94e-8 * erg2eV
    e_n2 = h*c / 5889.97e-8 * erg2eV
    u = part_func_na(temp)
    chi = [0, e_n1, e_n2]
    g = [2, 2, 4]
    relnrs = g[s]/u[r]*np.exp(-(chi[s])/(keV*temp))
    return relnrs


# Saha
def saha_na(temp, el_dens, ion_stage, e_ionization):
    kevt = keV*temp
    kergt = k_erg * temp
    u = part_func_na(temp)
    u = np.append(u, 2)  # append element to array
    saha_const = (2.*np.pi*m_e*kergt/(h*h)) ** (3./2) * 2. / el_dens
    n_stage = np.zeros(4)
    n_stage[0] = 1.
    for r in range(3):
        n_stage[r+1] = n_stage[r]*saha_const*u[r+1]/u[r] * np.exp(-e_ionization[r] / kevt)
    n_total = np.sum(n_stage)
    n_stage_rel = n_stage/n_total
    return n_stage_rel[ion_stage - 1]


# Saha Boltzmann
def sahabolt_na(temp, eldens, ionstage, level, e_ionization):
    return saha_na(temp, eldens, ionstage, e_ionization) * boltzmann_na(temp, ionstage, level)


# Doppler width testing
def doppler_width(wav, temp, v_t, m):
    # Takes in central wavelength in cm,temperature in K, v_t in km/s, and m in grams and returns dopplerwidth in cm.
    return wav/c*np.sqrt(2. * k_erg * temp / m + v_t * v_t * 1e10)


# Van der Waal broadening Na I D1 line
def rsq_nad(s, wav, e_ionization):  # put constants in statement here
    e_n = np.zeros(3, dtype=float)
    e_n[1] = h*c / wav[0] * erg2eV
    e_n[2] = h*c / wav[1] * erg2eV
    z = 1.
    rydberg = 13.6
    el = [0., 1., 1.]
    nstar_sq = rydberg * z**2 / (e_ionization[0] - e_n[s-1])
    rsq = nstar_sq / 2. / z**2 * (5*nstar_sq + 1 - 3*el[s-1]*(el[s-1] + 1))
    return rsq


def gammavdw_nad(temperature, p_gas, s, wav, e_ionization):
    rsq_u = rsq_nad(s, wav, e_ionization)
    rsq_l = rsq_nad(1, wav, e_ionization)
    log_vdw = 6.33 + 0.4 * np.log10(rsq_u - rsq_l) + np.log10(p_gas) - 0.7 * np.log10(temperature)
    return 10**log_vdw


# Define Voigt function
def voigt(a, u):
    # Calculates the Voigt function for values u and a
    z = u + 1.j*a
    return special.wofz(z).real
