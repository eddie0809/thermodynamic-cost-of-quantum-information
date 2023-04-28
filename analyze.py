#!/usr/bin/env python3

import glob
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colormaps as cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# for using tex formatting and font in plots

plt.rcParams.update({"text.usetex": True,}) 
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}\inputencoding{utf8}\usepackage{amsmath}\usepackage{amssymb}\usepackage{dsfont}\usepackage{mathtools}\usepackage{physics}'
mpl.rcParams['font.family'] = ['sans-serif']

cmap = cm['viridis']

files = [file for file in glob.glob("pickled_data/N=6/arr*")]
files.sort()

quantities = []
time_evo = []
integrated = []
discord = []
alpha = []

with open(files[0], 'rb') as f:
    the_dict = pickle.load(f)
    discord = the_dict[1]
    arr_time = the_dict[0]
    #discord.append(the_dict['discord'])
    #alpha.append(the_dict['alpha'])
    #quantities.append(the_dict['quantities'])
    #time_evo.append(the_dict['time_evo'])
    #integrated.append(the_dict['integrated'])

print(discord)
#t = [i/1000 for i in range(int(np.pi*1000)+1)]
fig, ax = plt.subplots(figsize=(14.5/2.54,6/2.54))
ax.plot(discord, np.array(arr_time)*1e-3, color=cmap(.5), label =r'Arrival Time')
#for ind, traj in enumerate(integrated[0]):
#    ax.plot(t, traj[:3142], label=r'$\expval{\sigma^z_{'+str(ind+1)+r'}}$', color=cmap((ind+1)/7))
#ax.set_xlim(-.05,np.pi+.05)

ax.set_xlabel(r'$D(\alpha)$')
ax.set_ylabel(r'Arrival time')
#
fig.tight_layout()
#fig.subplots_adjust(right=.83)
#fig.legend(loc='right')
plt.savefig('Review-2023-03-22/arr_time_over_discord.pdf')
plt.show()
#