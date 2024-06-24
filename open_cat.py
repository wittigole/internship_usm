import numpy as np

def mass_to_radius(m, delta = 500):
    #return the radius in Mpc given the mass in 10^14 Msun
    H_0 = 70*1e3/(3.086e22) #s^-1
    G = 6.67e-11 #m^3/kg/s^2
    rho_c = 3 * H_0**2 / (8 * np.pi * G)
    rho = delta * rho_c
    # m in kg, rho in kg/m^3
    r = (3 * m * 2e30 / (4 * np.pi * rho)) ** (1./3)
    return r #in m

#------------------------------------------------------------------

data = np.load('./new_mass_inon.npy')
m500_new = np.zeros(data.shape[1])

for i in range(data.shape[1]):
    m500_new[i] = float(data[1][i])
    """ print('Cluster: ' + str(data[0][i]) + ', M_halo: ' + str(float(data[1][i]) *1e-14) +\
           ' +/- ' + str(float(data[2][i]) * 1e-14) + ' 10^14 Msun') """
       
r500_new = mass_to_radius(m500_new, delta = 500)
c200 = 3.59
m200 = 6e14
r200 = mass_to_radius(m200, delta = 200)
c500 = np.median(r500_new) / r200 * c200
print(r500_new)
print('median r500: ' + str(np.median(r500_new)))
print('median r200: ' + str(r200))
print('median c200: ' + str(c200))
print('median c500: ' + str(c500))
