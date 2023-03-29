#!/usr/bin/env python3

import glob
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

plt.rcParams.update({"text.usetex": True,}) 
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}\inputencoding{utf8}\usepackage{amsmath}\usepackage{amssymb}\usepackage{dsfont}\usepackage{mathtools}\usepackage{physics}'
mpl.rcParams['font.family'] = ['serif']

cmap = mpl.cm.get_cmap('viridis')
cmap_2 = mpl.cm.get_cmap('plasma')


#npzfile = np.load("fid_of_discord.npz")
#npzfile_both = np.load("new_two_fid_of_discord.npz")
#
#fid_of_t_local = npzfile['fid_of_t']
#discord_local = npzfile['discord']
#
#fid_of_t_both_qubits = npzfile_both['fid_of_t']
#discord_both_qubits = npzfile_both['discord']
files = [file for file in glob.glob("pickled_data/N=6/new_two_*")]
files.sort()
files_2 = [file for file in glob.glob("pickled_data/N=6/metadata_for_*")]
files_2.sort()

fidelity = []
discord = []
t=[]
for ind, file in enumerate(files):
    with open(file, 'rb') as f:
        #local_dict = pickle.load(f)
        #print(type(local_dict))
        fidelity.append(pickle.load(f))
for file in files_2:
    with open(file, 'rb') as f:
        local_dict = pickle.load(f)
        discord.append(local_dict['discord'])
        t.append(local_dict['t'])


fid_max = np.zeros_like(discord)
t_of_fid_max = np.zeros_like(discord)

#for ind, fid in enumerate(fidelity.T):
#    fid_max[ind] = np.max(fid)
#    t_of_fid_max[ind] = np.argmax(fid)

#print(t_of_fid_max)

t_stop = np.pi
dt = 1e-3
t = np.arange(0, t_stop, dt)

fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
#for ind, discord_value in enumerate(discord_local):
#    ax.plot(t, fid_of_t_local[:,ind], '.', markersize=.5, zorder=ind, alpha=.6, zs=discord_value, zdir='y', color=cmap(discord_value/max(discord_local)))
    
for ind, discord_value in enumerate(discord):
    zs = discord_value
    ax.plot(t, fidelity[ind], '.', markersize=.5, zorder=ind, zs=discord_value, zdir='y', color=cmap_2(discord_value/max(discord)))
    #ax.scatter(t_of_fid_max[ind], fid_max[ind], zs, '.', color='black')

ax.set_xlim(0,3.2)
ax.set_ylim(0,0.125)
ax.set_zlim(1/2,1)
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$D(\alpha)$')
ax.set_zlabel(r'Fidelity')
ax.view_init(elev=0., azim=270)
plt.show()