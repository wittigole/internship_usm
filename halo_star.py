from __future__ import division
import pymultinest # type: ignore
import numpy as np # type: ignore
#from scipy.special import kn # type: ignore
#from scipy.special import iv # type: ignore
#import scipy.integrate as integrate # type: ignore
from astropy.table import Table # type: ignore
# from scipy.interpolate import interp1d # type: ignore
import astropy as ap # type: ignore
# from scipy import stats # type: ignore
from astropy.constants import G, c as splgt # type: ignore
from astropy.cosmology import FlatLambdaCDM, wCDM # type: ignore
from astropy import units as u # type: ignore
import json
import scipy # type: ignore
import os
import math
import h5py # type: ignore
from scipy.interpolate import UnivariateSpline # type: ignore
from astropy.cosmology import WMAP9 as cosmo1 # type: ignore
from scipy.interpolate import griddata # type: ignore
from scipy.special import gamma # type: ignore
import sys
import time

##### data #####

m500, m_star, z, m500_err, m_star_err= np.load('/home/o/Ole.Wittig/files/stellar_masses.npy') ### load your numpy file M500, M_star, redshift, M500_err, M_star_err
# data, data1, data2, data3, data4= np.load('/home/o/Ole.Wittig/files/stellar_masses.npy')

print(1)

mask = np.where(np.logical_and(m500 != 0, m_star != 0))[0]

data = m500[mask]
data1 = m_star[mask]
data2 = z[mask]
data3 = m500_err[mask]
data4 = m_star_err[mask]

# print(2)

# data, data1, data2, data3, data4 = np.load('/project/ls-mohr/users/aditya.singh/Ole_scaling_relation/halo_icm_redshift_combined_300_test_adi_all_red.npy')

# data = data[:60]
# data1 = data1[:60]
# data2 = data2[:60]
# data3 = data3[:60]
# data4 = data4[:60]

binval = 2 ### used for easily creating directories for saving the output, 1 is fine for the first run

#halo scaling relation usually used for fitting power laws (here in log-form)
#-> original form: M_star = A1*(M500/M_piv)**B1*((1+z)/(1+z_piv))**C1
def halo_star_scaling(M500,z,A1,B1,C1):
    # both pivot values are often median or mean values of the data
    M_piv = 4.8e14 #M_sun
    z_piv = 0.6 #redshift
    out = np.log(A1) + B1*np.log(M500/M_piv) + C1*np.log((1+z)/(1+z_piv))
    return out

#log-normal distribution: if X is normally distributed then Y = exp(X) is log-normal dostributed
# -> in this case m_star has a log-normal distribution and thus ln(m_star) is normally distributed
def lognormal(x1,lnM_star,D1):
    
    fac = 1/(np.exp(x1)*D1*np.sqrt(2*np.pi))
    fac1 = np.exp( -0.5*((x1-lnM_star)/D1)**2 )
    return np.log(fac*fac1)



class Priors(object):
    
  def __init__(self):
    pass
        
        
  def UniformPrior(self,r,x1,x2):
    """Uniform[0:1]  ->  Uniform[x1:x2]"""
    return x1+r*(x2-x1)
    
  def GaussianPrior(self,r,mu,sigma):
    """Uniform[0:1]  ->  Gaussian[mean=mu,variance=sigma**2]"""
    from math import sqrt
    from scipy.special import erfcinv  # type: ignore
    if (r <= 1.0e-16 or (1.0-r) <= 1.0e-16):
        return -1.0e32
    else:
        return mu+sigma*sqrt(2.0)*erfcinv(2.0*(1.0-r))
        
  def LogPrior(self,r,x1,x2):
    """Uniform[0:1]  ->  LogUniform[x1:x2]"""
    if (r <= 0.0):
            return -1.0e32
    else:
        from math import log10
        lx1=log10(x1); lx2=log10(x2)
        return 10.0**(lx1+r*(lx2-lx1))

def prior(cube, ndim, nparams):
    pri=Priors()
    
    ################ Uniform Priors ################
    
    #these are the starting values the MCMC algorithm will use for the parameters A1, B1, C1 and the scatter D1

    cube[0] = pri.UniformPrior(cube[0],1e12,2e13)  #Adjust the range accroding to you M_star, should be wide enough
    cube[1] = pri.UniformPrior(cube[1],0.1,1.5) #exponent of the power law, around 0.8 for my sample
    cube[2] = pri.UniformPrior(cube[2],-3,3) #redshift range, around 0.5 for my sample
    cube[3] = pri.UniformPrior(cube[3],0.001,1) #scatter

    
def loglikelihood(cube, ndim, nparams):
 
    A1 = cube[0]
    B1 = cube[1]
    C1 = cube[2]
    D1 = cube[3]
        
    LL = 0
    for i in range(len(data)):
        data_M500 = data[i]
        data_lnM_star = np.log(data1[i])
        z = data2[i]
        ln_m500_err = data3[i]/data[i]   ### measurement scatter in ln_m500
        ln_mstar_err = data4[i]/data1[i]  ### measurement scatter in ln_m_star

        D_tot = np.sqrt(D1**2 + ln_mstar_err**2 + B1**2*ln_m500_err**2)    #### Total scatter

        pred_lnM_star = halo_star_scaling(data_M500, z, A1, B1, C1) 
        LL += lognormal(data_lnM_star, pred_lnM_star, D_tot)
    LL_final = LL


    return LL_final

print(3)

parameters = ["A1","B1","C1","D1"]

n_params = len(parameters)

directory = '/home/o/Ole.Wittig/files/chain_no_measurement_error_%s'%binval  #### change it to your directory
if not os.path.exists(directory):
    os.mkdir(directory)

datafile = '/home/o/Ole.Wittig/files/chain_no_measurement_error_%s/'%binval  #### change it to your directory

##### Multinest starts here, no need to change the parameters!
print(4)
start = time.time()

pymultinest.run(loglikelihood, prior, n_params, outputfiles_basename=datafile + '_1_', resume = False, verbose = True, n_live_points= 400, evidence_tolerance = 0.3, init_MPI = False, sampling_efficiency= 0.3,const_efficiency_mode = False, log_zero=-1e+90)
json.dump(parameters, open(datafile + '_1_params.json', 'w'))

a = pymultinest.Analyzer(outputfiles_basename=datafile + '_1_', n_params = n_params)

end = time.time()

print('MCMC curve fit took ', np.round(end-start,1), 'seconds.')

