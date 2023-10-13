import pickle
import numpy as np
import datetime
import time
from qutip import *
import matplotlib.pyplot as plt
from matplotlib import colormaps as cm

import xy_chain

plt.rcParams['text.usetex'] = True  # to use LaTeX in figures
plt.rcParams["text.latex.preamble"]= r'\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern, mathpazo}\inputencoding{utf8}\usepackage{amsmath}\usepackage{amssymb}\usepackage{dsfont}\usepackage{mathtools}\usepackage{physics}\usepackage{siunitx}'
plt.rcParams['font.family'] = ['serif']

cmap = cm['magma']

beta = [1e-3, 1e-3]
N = 5
_lambda = 2*np.sqrt(N)/N
t_stop = np.pi/_lambda
dt = 5e-4
t = np.arange(0, t_stop, dt)

system_0 = xy_chain.HeisenbergXY(t, N, 0, beta, corr='therm', lamb=_lambda)
system_corr = xy_chain.HeisenbergXY(t, N, 1, beta, corr='therm', lamb=_lambda)

system_0.e_dot_test()
system_corr.e_dot_test()

etot = system_corr.quantities['e_dot_test']
edot_0 = system_0.quantities['e_dot_test']

edot_chi = etot-edot_0

plt.plot(edot_chi)
plt.plot(edot_0)
plt.plot(etot)
edotdot_chi = np.diff([0] + edot_chi)/dt
print(edotdot_chi)
plt.plot(edotdot_chi)
plt.plot([0 for _ in etot], color='grey')
ind = np.argmin(edotdot_chi)
print(ind, edot_chi[ind], etot[ind])
plt.show()
