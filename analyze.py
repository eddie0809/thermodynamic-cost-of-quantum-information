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

files = [file for file in glob.glob("pickled_data/N=6/dict_*")]
files.sort()

quantities = []
time_evo = []
integrated = []
discord = []
alpha = []
for ind, file in enumerate(files):
    with open(file, 'rb') as f:
        the_dict = pickle.load(f)
        discord.append(the_dict['discord'])
        alpha.append(the_dict['alpha'])
        quantities.append(the_dict['quantities'])
        time_evo.append(the_dict['time_evo'])
        integrated.append(the_dict['integrated'])

arr_time = [1e-3*quantities[ind]['local minima in rel ent'][0][0] for ind in range(len(files))]
#print(discord, arr_time)
fig, ax = plt.subplots(figsize=(12,5))
ax.plot(discord, arr_time, '.')
ax.set_xscale('log')
#ax.plot(discord,alpha)
plt.show()
