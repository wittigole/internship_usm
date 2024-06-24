import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
import funcs
import h5py # type: ignore

f = h5py.File('./stellar_masses.hdf5', 'r')
z = f['z'][:]
m500 = f['m500'][:]
m500_err = f['m500_err'][:]
m500_old = f['m500_old'][:] * 1e14
m500_old_err = f['m500_old_err'][:] * 1e14
m_star = f['m_star'][:] * 1e12
m_star_old = f['m_star_old'][:] * 1e12

f.close()

mask = np.where(m500 != 0)[0]

low_z = np.where(z < 0.3)[0]
mid_z = np.where(np.logical_and(z > 0.3, z < 0.6))[0]
high_z = np.where(z > 0.6)[0]

# fig, ax = plt.subplots(1,1, figsize=(8,6))
# ax.errorbar(m500[low_z], m_star[low_z], xerr=m500_err[low_z], fmt='o', color = 'gold', label='z < 0.3')
# ax.errorbar(m500[mid_z], m_star[mid_z], xerr=m500_err[mid_z], fmt='o', color = 'darkred', label='0.3 < z < 0.6')
# ax.errorbar(m500[high_z], m_star[high_z], xerr=m500_err[high_z], fmt='o', color = 'black', label='z > 0.6')
# ax.legend(loc = 'lower right')
# ax.set(xscale='linear', yscale='linear', xlabel='M500 [Msun]', ylabel='Mstar [Msun]')
#fig.set_tight_layout(True)
# plt.show()
#fig.savefig('./stellar_masses.eps')


# fig, ax = plt.subplots(1,3, figsize=(18,6))
# ax[0].scatter(m500[low_z], m_star[low_z], color = 'gold', label='new', marker = 'o')
# ax[0].scatter(m500_old[low_z], m_star_old[low_z], color = 'gold', label='old', marker = 'd')

# for i in range(low_z.shape[0]):
#     ax[0].plot([m500[low_z[i]], m500_old[low_z[i]]], [m_star[low_z[i]], m_star_old[low_z[i]]], color = 'gold',\
#                 alpha = 0.3)

# ax[0].title.set_text('z < 0.3')

# ax[1].scatter(m500[mid_z], m_star[mid_z], color = 'darkred', label='new', marker = 'o')
# ax[1].scatter(m500_old[mid_z], m_star_old[mid_z], color = 'darkred', label='old', marker = 'd')
# for i in range(mid_z.shape[0]):
#     ax[1].plot([m500[mid_z[i]], m500_old[mid_z[i]]], [m_star[mid_z[i]], m_star_old[mid_z[i]]], color = 'darkred',\
#                 alpha = 0.3)
# ax[1].title.set_text('0.3 < z < 0.6')

# ax[2].scatter(m500[high_z], m_star[high_z], color = 'black', label='new', marker = 'o')
# ax[2].scatter(m500_old[high_z], m_star_old[high_z], color = 'black', label='old', marker = 'd')
# for i in range(high_z.shape[0]):
#     ax[2].plot([m500[high_z[i]], m500_old[high_z[i]]], [m_star[high_z[i]], m_star_old[high_z[i]]], color = 'black',\
#                 alpha = 0.3)
# ax[2].title.set_text('z > 0.6')


# ax[2].legend(loc = 'lower right')
# for i in range(3):
    
#     ax[i].set(xscale='linear', yscale='linear', xlabel='M500 [Msun]', ylabel='Mstar [Msun]', ylim = (0,1e13))

# fig.set_tight_layout(True)
# plt.show()
#fig.savefig('./stellar_masses_comp.eps')

# fig, ax = plt.subplots(1,2, figsize=(12,6))
# ax[0].scatter(m500[mask], m_star[mask]/m_star_old[mask], color = 'black', marker = 'o')
# ax[0].set(xscale='log', xlabel='M500 [Msun]', ylabel='Mstar_new/Mstar_old', xlim = (3e14,1.2e15))
# ax[0].hlines(1, 3e14, 1.2e15, color='gray', linestyle='--', linewidth = 2)

# ax[1].scatter(z[mask], m_star[mask]/m_star_old[mask], color = 'black', marker = 'o')
# ax[1].set(xlabel='z', ylabel='Mstar_new/Mstar_old', xlim = (0.15,1.05))
# ax[1].hlines(1, 0.1, 1.1, color='gray', linestyle='--', linewidth = 2)

# fig.set_tight_layout(True)
# plt.show()
#fig.savefig('./stellar_mass_ratio_distribution.eps')

# fig, ax = plt.subplots(1,1, figsize=(8,6))
# ax.scatter(m500[mask]/m500_old[mask], m_star[mask]/m_star_old[mask], color = 'black', marker = 'o')
# ax.set(xlabel='M500_new/M500_old', ylabel='Mstar_new/Mstar_old', xlim = (0.95,1.6), ylim = (0.99,1.2))
# x = np.linspace(0.9,1.7,200)
# ax.plot(x, x, color = 'gray', linestyle='--', linewidth = 2)

# fig.set_tight_layout(True)
# plt.show()
#fig.savefig('./stellar_vs_halo_mass_ratio.eps')

fig, ax = plt.subplots(1,1, figsize=(8,6))
ax.errorbar(m500[mask], m_star[mask], xerr=m500_err[mask], fmt='o', color = 'gold', label='data')

chiu = np.argsort(m500[mask])
m500_chiu = m500[mask[chiu]]
z_chiu = z[mask[chiu]]
ax.plot(m500_chiu, funcs.obs_mass_z_relation(m500 = m500_chiu, z = z_chiu, m_piv = 4.8e14, z_piv = 0.6, a = 4e12, b = 0.8, c = 0.05), color = 'black', label='Chiu+ 2018', linewidth = 2)

m500_tng = np.array([2.34826599e+14, 3.68474438e+14, 4.62433626e+14, 6.15112331e+14, 7.03479337e+14])
m_star_tng = np.array([4.90334742e+12, 7.30638582e+12, 8.26568553e+12, 1.19388357e+13, 1.29396071e+13])
ax.plot(m500_tng, m_star_tng, color = 'red', label='TNG300 (incl. ICL)', linewidth = 2)

ax.legend(loc = 'lower right')
ax.set(xscale='log', yscale='log', xlabel='M500 [Msun]', ylabel='Mstar [Msun]', xlim = (3e14,1.5e15), ylim = (1e12,1.5e13))

fig.set_tight_layout(True)
plt.show()
fig.savefig('./plots/stellar_vs_halo_mass.jpg', format = 'jpg')