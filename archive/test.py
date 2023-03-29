#!/usr/bin/env python3

import glob
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


files_2 = [file for file in glob.glob("pickled_data/N=6/2022-12-16/quantities*")]
files_2.sort()

discord = []
rho_T = []

for ind, file in enumerate(files_2):
    with open(file, 'rb') as f:
        local_dict = pickle.load(f)
#        #i_dot_sq_max.append(local_dict['i_dot_sq_max'])
        discord.append(local_dict['discord'])
        rho_T.append()

