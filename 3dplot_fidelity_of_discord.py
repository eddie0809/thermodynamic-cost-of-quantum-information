#!/usr/bin/env python3

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

plt.rcParams.update({"text.usetex": True,}) 
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}\inputencoding{utf8}\usepackage{amsmath}\usepackage{amssymb}\usepackage{dsfont}\usepackage{mathtools}\usepackage{physics}'
mpl.rcParams['font.family'] = ['serif']

cmap = mpl.cm.get_cmap('viridis')
cmap_2 = mpl.cm.get_cmap('plasma')


npzfile = np.load("fid_of_discord.npz")
npzfile_both = np.load("two_fid_of_discord.npz")

fid_of_t_local = npzfile['fid_of_t']
discord_local = npzfile['discord']

fid_of_t_both_qubits = npzfile_both['fid_of_t']
discord_both_qubits = npzfile_both['discord']

t_stop = np.pi
dt = 1e-3
t = np.arange(0, t_stop, dt)

fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
for ind, discord_value in enumerate(discord_local):
    ax.plot(t, fid_of_t_local[:,ind], '.', markersize=.5, zorder=ind, alpha=.6, zs=discord_value, zdir='y', color=cmap(discord_value/max(discord_local)))
    
for ind, discord_value in enumerate(discord_both_qubits):
    ax.plot(t, fid_of_t_both_qubits[:,ind], '.', markersize=.5, zorder=ind, zs=discord_value, zdir='y', color=cmap_2(discord_value/max(discord_both_qubits)))
ax.set_xlim(0,3.2)
ax.set_zlim(1/2,1)
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$D(\alpha)$')
ax.set_zlabel(r'Fidelity')
ax.view_init(elev=0., azim=270)
plt.show()