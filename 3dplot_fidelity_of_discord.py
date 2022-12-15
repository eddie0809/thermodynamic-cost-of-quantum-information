#!/usr/bin/env python3

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

plt.rcParams.update({"text.usetex": True,}) 
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}\inputencoding{utf8}\usepackage{amsmath}\usepackage{amssymb}\usepackage{dsfont}\usepackage{mathtools}\usepackage{physics}'
mpl.rcParams['font.family'] = ['serif']

cmap = mpl.cm.get_cmap('viridis')

npzfile = np.load("fid_of_discord.npz")

fid_of_t = npzfile['fid_of_t']
discord = npzfile['discord']

t_stop = np.pi
dt = 1e-3
t = np.arange(0, t_stop, dt)

fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
for ind, discord_value in enumerate(discord):
    ax.plot(t, fid_of_t[:,ind], '.', markersize=.5, zs=discord_value, zdir='y', color=cmap(discord_value/max(discord)))
ax.set_xlim(0,3.2)
ax.set_zlim(1/np.sqrt(2),1)
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$D(\alpha)$')
ax.set_zlabel(r'Fidelity')
ax.view_init(elev=25., azim=-90)
plt.show()