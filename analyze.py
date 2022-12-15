#!/usr/bin/env python3

import os
import glob
import pickle
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.optimize import curve_fit

# for using tex formatting and font in plots

plt.rcParams.update({"text.usetex": True,}) 
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}\inputencoding{utf8}\usepackage{amsmath}\usepackage{amssymb}\usepackage{dsfont}\usepackage{mathtools}\usepackage{physics}'
mpl.rcParams['font.family'] = ['serif']

cmap = mpl.cm.get_cmap('viridis')

files = [file for file in glob.glob("pickled_data/N=6/2022-12-15/quantities/*")]
files.sort()
i_dot_sq_max = []
discord = []
fid_of_t = []
for file in files:
    with open(file, 'rb') as f:
        local_dict = pickle.load(f)
        i_dot_sq_max.append(local_dict['i_dot_sq_max'])
        discord.append(local_dict['discord'])
        fid_of_t.append(local_dict['fidelity'])


plt.plot(discord, i_dot_sq_max/i_dot_sq_max[-1], '.', color=cmap(1/2), label = r'$\frac{\dot{I}_\mathrm{max}^2(D(\alpha))}{\dot{I}_\mathrm{max}^2(D(\alpha=0))}$')

#t_stop = np.pi
#dt = 1e-3
#t = np.arange(0, t_stop, dt)
##
#fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
#for ind, discord_value in enumerate(discord):
#    ax.plot(t, fid_of_t[ind], '.', markersize=.5, zs=discord_value, zdir='y', color=cmap(discord_value/max(discord)))
#ax.set_xlim(0,3.2)
#ax.set_zlim(1/np.sqrt(2),1)
#ax.view_init(elev=20., azim=-90)
plt.show()

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