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

files = [file for file in glob.glob("pickled_data/N=6/new_two_*")]
files.sort()
files_2 = [file for file in glob.glob("pickled_data/N=6/metadata_for_*")]
files_2.sort()

fidelity = []
discord = []
t=[]

for ind, file in enumerate(files):
    with open(file, 'rb') as f:
        fidelity.append(pickle.load(f))
for file in files_2:
    with open(file, 'rb') as f:
        local_dict = pickle.load(f)
        discord.append(local_dict['discord'])
        t.append(local_dict['t'])