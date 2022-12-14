#!/usr/bin/env python3

import glob
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# for using tex formatting and font in plots

plt.rcParams.update({"text.usetex": True,}) 
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}\inputencoding{utf8}\usepackage{amsmath}\usepackage{amssymb}\usepackage{dsfont}\usepackage{mathtools}\usepackage{physics}'
mpl.rcParams['font.family'] = ['serif']

cmap = mpl.cm.get_cmap('viridis')

files = [file for file in glob.glob("pickled_data/N=6/two_*")]
files.sort()
files_2 = [file for file in glob.glob("pickled_data/N=6/quantities*")]
files_2.sort()
#i_dot_sq_max = []
discord = []

for ind, file in enumerate(files_2):
    with open(file, 'rb') as f:
        local_dict = pickle.load(f)
        #i_dot_sq_max.append(local_dict['i_dot_sq_max'])
        discord.append(local_dict['discord'])

fid_of_t = np.zeros((3142,len(files)))
i_dot_sq_of_t = np.zeros((3141, len(files)))
for ind, file in enumerate(files):
    with open(file, 'rb') as f:
        local_dict = pickle.load(f)
        #i_dot_sq_max.append(local_dict['i_dot_sq_max'])
        #discord.append(local_dict['discord'])
        fid_of_t[:,ind] = local_dict
        #i_dot_sq_of_t[:,ind] = local_dict['i_dot_sq']
#discord = np.array(discord)
print(fid_of_t)
with open("two_fid_of_discord.npz", 'wb') as f:
    np.savez(f, fid_of_t=fid_of_t, discord=discord)
#
#fig, ax =plt.subplots(figsize=(8, 4.5))
#ax.plot(discord, i_dot_sq_max/i_dot_sq_max[-1], '.', color=cmap(1/2), label = r'$\frac{\dot{I}_\mathrm{max}^2(D(\alpha))}{\dot{I}_\mathrm{max}^2(D(\alpha=0))}$')
#t_stop = np.pi
#dt = 1e-3
#t = np.arange(0, t_stop, dt)
###
#fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
#for ind, discord_value in enumerate(discord):
#    ax.plot(t[:-1], i_dot_sq_of_t[:,ind], '.', markersize=.5, zs=discord_value, zdir='y', color=cmap(discord_value/max(discord)))
#    new_zs = np.argmax(i_dot_sq_of_t[:,ind])*1e-3
#    print(new_zs)
#    ax.scatter(new_zs, i_dot_sq_max[ind], zdir='y',zs=discord_value, color='black')
#ax.set_ylim(0,max(discord))
#ax.set_ylabel(r"$D(\alpha)$")
#ax.set_xlabel(r"$t$")
##ax.set_ylabel(r"$\frac{\dot{I}_\mathrm{max}^2(D(\alpha))}{\dot{I}_\mathrm{max}^2(D(\alpha=0))}$")
##ax.set_zlim(1/np.sqrt(2),1)
#ax.view_init(elev=36., azim=70)
#plt.savefig("plots/i_dot_sq_of_discord_3dplot_with_max.pdf")
#plt.show()

#fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12,8))

#for n in range(N+1):
#    ax.plot(t, sz[n], label=r"$\expval{\sigma^z_{%d}}$"%(n+1), lw=1.5, color=cmap((n+1)/(N+2)))
#ax.set_xlabel(r'Time')
#ax.set_ylabel(r'$\expval{\sigma^z_i}$')
#ax1.plot(t, fidelity_of_t, label=r"$F(\rho_\mathrm{corr, 12}, \rho_{56}^t)$", color=cmap(1/3))
#ax2.plot(t, rel_ent_of_t, label=r'$\mathcal{S}(\rho_\mathrm{corr, 12} \mid\mid \rho_{56}^t)$', color=cmap(2/3))
#ax2.set_yscale('log')
#ax.axhline(0, color='grey', ls='--')
#ax1.axhline(1, color='grey', ls='--')
#ax2.plot(t[:-1], i_dot_sq_perf, color=cmap(2/3))
#ax1.legend()
#fig.legend(loc=8, ncol=4)
#fig.tight_layout()
#fig.subplots_adjust(bottom=.1)
#plt.show()