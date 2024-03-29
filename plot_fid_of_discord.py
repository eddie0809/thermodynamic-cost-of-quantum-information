#!/usr/bin/env python3

import glob
import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib import colormaps as cm
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.mplot3d import Axes3D
from qutip import *


plt.rcParams.update({"text.usetex": True, })
plt.rcParams['text.latex.preamble'] = r"""\usepackage[utf8]{inputenc}
    \usepackage[T1]{fontenc}\usepackage{lmodern}\inputencoding{utf8}
    \usepackage{amsmath}\usepackage{amssymb}\usepackage{dsfont}
    \usepackage{mathtools}\usepackage{physics}
    """
mpl.rcParams['font.family'] = ['sans-serif']

cmap = cm['RdBu']
cmap_2 = cm['plasma']


files = [file for file in glob.glob("pickled_data/N=6/new_two_*")]
files.sort()
files_2 = [file for file in glob.glob("pickled_data/N=6/metadata_for_*")]
files_2.sort()

integration = [
    file for file in glob.glob("pickled_data/N=6/system_integrated*")
    ]
integration.sort()

N = 5
two_fidelity = []
integrated = []
discord = []
t = []

for ind, file in enumerate(files):
    with open(file, 'rb') as f:
        two_fidelity.append(pickle.load(f))

for file in files_2:
    with open(file, 'rb') as f:
        local_dict = pickle.load(f)
        discord.append(local_dict['discord'])
        t.append(local_dict['t'])

for file in (integration[0], integration[-1]):
    with open(file, 'rb') as f:
        integrated.append(pickle.load(f))

t = t[0]
time_evo_zero_alpha = integrated[0]
time_evo_max_alpha = integrated[1]

one_fidelity_0 = [
    fidelity(time_evo_zero_alpha[0].ptrace(0), time_evolved.ptrace(N))
    for time_evolved in time_evo_zero_alpha
    ]
one_fidelity_1 = [
    fidelity(time_evo_max_alpha[0].ptrace(0), time_evolved.ptrace(N))
    for time_evolved in time_evo_max_alpha
    ]


fig, ((ax00, ax01), (ax10, ax11)) = plt.subplots(
    nrows=2, ncols=2, sharex=True, figsize=(14.5/2.54, 14.5/2.54)
    )
ax00.plot(t, two_fidelity[0],
          color=cmap(3/4),
          label=r'$D(\alpha) = {disc}$'.format(disc=discord[0]))
ax00.set_ylabel(r"Fidelity")
ax01.plot(t, two_fidelity[-1],
          color=cmap(1/4),
          label=r'$D(\alpha) = {disc}$'.format(disc=discord[-1]))
ax10.plot(t, one_fidelity_0,
          color=cmap(3/4),
          label=r'$D(\alpha) = {disc}$'.format(disc=discord[0]))
ax10.set_xlabel(r'$t$')
ax10.set_ylabel(r"Fidelity")
ax11.plot(t, one_fidelity_1,
          color=cmap(1/4),
          label=r'$D(\alpha) = {disc}$'.format(disc=discord[-1]))
ax11.set_xlabel(r'$t$')
at = AnchoredText(r"$F(\rho_{12}(0), \rho_{N,N-1}(t))$",
                  prop=dict(size=11), frameon=True, loc='lower right')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.05")
ax00.add_artist(at)
# ax01.add_artist(at)
at0 = AnchoredText(r"$F(\rho_{1}(0), \rho_{N}(t))$",
                   prop=dict(size=11), frameon=True, loc='lower right')
at0.patch.set_boxstyle("round,pad=0.,rounding_size=0.05")
ax10.add_artist(at0)
# ax11.add_artist(at0)
fig.subplots_adjust(bottom=.22)
fig.legend(loc='lower center', ncols=4)
fig.tight_layout()
plt.savefig("plots/max_und_min_discord_fidelities.pdf")
plt.show()
