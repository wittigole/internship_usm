import numpy as np # type: ignore
from scipy.integrate import quad # type: ignore
from scipy.optimize import root, root_scalar # type: ignore

#define global constants

# physical constants
G = 6.67e-11 #m^3/kg/s^2, gravitational constant
Msun = 1.989e30 #kg, solar mass
Mpc = 3.086e22 #m, mega parsec
c = 3e8 #m/s, speed of light

# cosmological parameters
H_0 = 70*1e3/Mpc #s^-1, Hubble constant
omega_m = 0.315 #matter density parameter
omega_l = 0.685 #dark energy density parameter

def dimless_hubble_param(z):
    return np.sqrt(omega_m * (1+z)**3 + omega_l)

def nfw_solution(r, c):
    """ return the integrated NFW profile up to r given the concentration. normalization is arbitrary and r is assumed to
     fulfill r = r_s * c """
    return (r/c)**3 * (np.log(1+c) - c/(1+c))

def nfw_solution_err(r, c, r_err, c_err):
    """ return the error on the integrated NFW profile up to r given the concentration. normalization is arbitrary and r is assumed to
     fulfill r = r_s * c """
    t1 = 3 * r**2 / c**3 * (np.log(1+c) - c/(1+c)) * r_err
    t2 = r**3 / c**4 * (4*c**2 + 3*c - 3*(1 + c)**2 * np.log(1+c)) / (1+c)**2 * c_err
    return np.sqrt(t1**2 + t2**2)

def mass_to_radius(m, m_err = None, delta = 500):
    """ return the radius in Mpc given the mass in Msun """
    rho_c = 3 * H_0**2 / (8 * np.pi * G)
    # rho_c *= dimless_hubble_param(z)**2
    rho = delta * rho_c
    # m in kg, rho in kg/m^3
    r = (3 * m * Msun / (4 * np.pi * rho)) ** (1./3)
    r_err = (3 * Msun/ (4 * np.pi * rho)) ** (1./3) * (1/3) * m**(-2/3) * m_err
    return r, r_err #in m

def mass_from_profile(r, c, r_err, c_err):
    """ return the mass of the galaxy cluster within r500 or r200 or whatever """
    rho_c = 3 * H_0**2 / (8 * np.pi * G)
    # rho_c *= dimless_huubble_param(z)
    rho = 500 * rho_c
    M = 4 * np.pi * rho * nfw_solution(r,c) #integrated NFW profile (3D)
    M_err = 4 * np.pi * rho * nfw_solution_err(r, c, r_err, c_err)

    return M / Msun, M_err / Msun #in solar masses

def value_from_radii(r200_guess, r500, c):
    """ return THE value at r200 (guess) given r500 and the concentration """
    num = (r200_guess/r500)**3 * (np.log((r200_guess + c*r500)/r200_guess) - c*r500 / (r200_guess + c*r500))
    den = np.log(1+c) - c/(1+c)
    value = num/den - 2.5
    return value

def find_r200(r500, c):
    """ return the r200 given r500 and the concentration """
    if isinstance(r500, (int, float)):
        r500 = np.array([r500])

    r200 = np.zeros(r500.shape[0])
    r200_guess = r500

    for i in range(r500.shape[0]):
        r200[i] = root(value_from_radii, args=(r500[i], c), x0=r200_guess[i]).x
        # r200[i] = root_scalar(value_from_radii, args=(r500[i], c), x0=r200_guess[i]).root
    return r200

def update_ratio(r500_old, r500_new, c500_new, c500_old, r500_old_err, r500_new_err, c500_new_err, c500_old_err):
    """ return the ratio between the old and the new stellar mass """
    num = nfw_solution(r500_new, c500_new) #new
    den = nfw_solution(r500_old, c500_old) #old

    num_err = nfw_solution_err(r500_new, c500_new, r500_new_err, c500_new_err)
    den_err = nfw_solution_err(r500_old, c500_old, r500_old_err, c500_old_err)

    ratio = num / den
    ratio_err = ratio * np.sqrt((num_err/num)**2 + (den_err/den)**2)
    return ratio, ratio_err

def integrand(z):
    # this is actually the inverse of the dimensionless Hubble parameter E(z) = H(z)/H_0 = sqrt(omega_m * (1+z)^3 + omega_l)
    return 1 / np.sqrt(omega_m * (1+z)**3 + omega_l)

def angle_to_size(z, angle):
    """ return the size in m given the redshift and the angle in arcmin """

    if isinstance(z, (int, float)):
        z = np.array([z])

    x = np.zeros(z.shape[0])
    for i in range(z.shape[0]):
        x[i], _ = quad(integrand, 0, z[i], args=())
    D_A = c/H_0 / (1+z) * x
    size = D_A * angle * np.pi / 10800 #in m
    return size #in m

def obs_mass_z_relation(m500, z, a, b, c):
    """scaling relation between halo mass, redshift and observable"""
    m_piv = 4.8e14
    z_piv = 0.6
    return a*(m500/m_piv)**b * ((1+z)/(1+z_piv))**c
