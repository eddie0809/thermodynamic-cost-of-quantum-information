#!/usr/bin/env python3

import os
import glob
import pickle
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# for using tex formatting and font in plots

plt.rcParams.update({"text.usetex": True,}) 
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}\inputencoding{utf8}\usepackage{amsmath}\usepackage{amssymb}\usepackage{dsfont}\usepackage{mathtools}\usepackage{physics}'
mpl.rcParams['font.family'] = ['serif']

cmap = mpl.cm.get_cmap('viridis')

files = [file for file in glob.glob("pickled_data/N=6/test_*")]
files.sort()
i_dot_sq_max = []
discord = []
arr_time = []
for file in files:
    with open(file, 'rb') as f:
        local_dict = pickle.load(f)
        i_dot_sq_max.append(local_dict['i_dot_sq_max'])
        discord.append(local_dict['discord'])
        arr_time.append(np.nanmin(local_dict['time_where_kld_isclose_0']))
print(arr_time)

plt.plot(discord, i_dot_sq_max)
plt.plot(discord[0:4], arr_time[0:4])
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