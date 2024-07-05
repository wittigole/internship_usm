import numpy as np
import funcs
import h5py

def delete_minus(valstr):
    return float(valstr.decode("utf-8").replace('-','0'))

cols = (2,3,5,9,11)
converters = {i: delete_minus for i in cols}

z = np.loadtxt('tab1.txt', skiprows=1, usecols=(2,))
tab2 = np.loadtxt('tab2.txt', skiprows=1, usecols=cols, converters = converters)

# r500 = tab2[:,0] #in arcmin
# r500 = funcs.angle_to_size(z,r500) #in Mpc

m500_old = tab2[:,1] * 1e14 #in Msun
m500_old_err = tab2[:,2] * 1e14 #in Msun
m_star_old = tab2[:,3] * 1e12 #in Msun
m_star_old_err = tab2[:,4] * 1e12 #in Msun

r500_old, r500_old_err = funcs.mass_to_radius(m500_old, m500_old_err, delta = 500) #in m

data = np.load('./new_mass_inon.npy')
m500_new = np.zeros(data.shape[1])
m500_new_err = np.zeros(data.shape[1])

z_new = np.zeros(data.shape[1])

for i in range(data.shape[1]):
    m500_new[i] = float(data[1][i])
    m500_new_err[i] = float(data[2][i])
    z_new[i] = float(data[3][i])

r500_new, r500_new_err = funcs.mass_to_radius(m500_new, m500_new_err, delta = 500)
c = 3.59
c_err = 0.2

data200_1, data200_2, data200_3, data200_4 = np.load('./new_mass_inon_M200c.npy')
# print(data200_2)
# print(data200_3)
# print(data200_4)
m200_new = data200_2.astype(float)
m200_new_err = data200_3.astype(float)
z_alt = data200_4.astype(float)

assert np.all(z_alt == z_new), 'wrong redshifts'

r200_new, r200_new_err = funcs.mass_to_radius(m200_new, m200_new_err, delta = 200)

r200_new_alt = funcs.find_r200(r500_new, c)

print(np.abs(r200_new - r200_new_alt)/r200_new)

c500_new = c * r500_new / r200_new
c500_new_err = c500_new * np.sqrt((c_err/c)**2 + (r500_new_err/r500_new)**2 + (r200_new_err/r200_new)**2)
#r200_old = funcs.find_r200(r500_old, c)
c500_old = c * r500_old / r200_new
c500_old_err = c500_old * np.sqrt((c_err/c)**2 + (r500_old_err/r500_old)**2 + (r200_new_err/r200_new)**2)

ratio, ratio_err = funcs.update_ratio(r500_old, r500_new, c500_new, c500_old, r500_old_err, r500_new_err, c500_new_err, c500_old_err)

mask = np.where(np.isfinite(ratio))[0]

m_star_new = m_star_old * ratio
m_star_new_err = m_star_new * np.sqrt((m_star_old_err/m_star_old)**2 + (ratio_err/ratio)**2)

# print('mass ratio \t sqrt(r/r_old)')
# for i in range(len(ratio)):
    # print(str(ratio[i]) + ' \t ' + str(np.sqrt(r500_new[i])/np.sqrt(r500_old[i])))

f = h5py.File('./files/stellar_masses.hdf5', 'w')
f.create_dataset('z', data=z_new)
f.create_dataset('m500_new', data=m500_new)
f.create_dataset('m500_new_err', data=m500_new_err)
f.create_dataset('m500_old', data=m500_old)
f.create_dataset('m500_old_err', data=m500_old_err)
f.create_dataset('r500_new', data=r500_new)
f.create_dataset('r500_new_err', data=r500_new_err)
f.create_dataset('r500_old', data=r500_old)
f.create_dataset('r500_old_err', data=r500_old_err)

f.create_dataset('r200', data=r200_new)
f.create_dataset('c500_new', data=c500_new)
f.create_dataset('c500_new_err', data=c500_new_err)

f.create_dataset('m_star_new', data=m_star_new)
f.create_dataset('m_star_err_new', data=m_star_new_err)
f.create_dataset('m_star_old', data=m_star_old)
f.create_dataset('m_star_err_old', data=m_star_old_err)
f.close()


np.save('files/stellar_masses.npy', (m500_new, m_star_new, z_new, m500_new_err, m_star_new_err))

#print(m500_new[0])
#print(m500_new_err[0])

#print(m_star_new[0])
#print(m_star_new_err[0])

#print(np.sqrt(r500_new[mask])/np.sqrt(r500[mask]))
#print(r500_new[mask]/r500[mask])
#print(m500_new[mask]/m500[mask]/1e14)
#print(c500_new[mask]/c500_old[mask])
#print(ratio[mask])
