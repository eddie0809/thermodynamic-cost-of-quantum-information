import pickle
import numpy as np
import datetime
import time
from qutip import *
import matplotlib as mpl
mpl.use('pgf')
import matplotlib.pyplot as plt
from matplotlib import colormaps as cm

import xy_chain

plt.rcParams.update({
    "font.family":'sans-serif',
    "font.size":8,
    "text.usetex":True,
    "pgf.texsystem": "lualatex",
    "pgf.rcfonts":False,
    "pgf.preamble": "\n".join([
         r"\usepackage{fontspec}",
         r"\setmainfont{Latin Modern Roman}",
         r"\usepackage{amsmath}",
         r"\usepackage{amssymb}",
         r"\usepackage{mathtools}", 
         r"\usepackage{dsfont}",
         r"\usepackage{physics}"
    ]),
})

cmap = cm['magma']

beta = [1e-3, 1e-3]
N = 5
_lambda = 2*np.sqrt(N)/N
t_stop = np.pi/_lambda
dt = 5e-4
t = np.arange(0, t_stop, dt)
testing = time.time()
print(time.time()-testing)
x_state = Qobj([
    [.4, 0, 0, 1.j*.04],
    [0, .3, 1.j*.06, 0],
    [0, -1.j*.06, .2, 0],
    [-1.j*.04, 0, 0, .1]
])

poss_alpha = [(i/100)**2 for i in range(101)]

arr_times = []
discord = []
crit_rho = {}
for alph in poss_alpha:
    system = xy_chain.HeisenbergXY(t, N, alph, beta, corr='therm', lamb=_lambda)
    print(f"system: {time.time()-testing}")
    system.arrival_time()
    print(f"arrival time: {time.time()-testing}")
    ind_arr = np.min(system.quantities["local minima in rel ent"])
    the_arr_time = t[ind_arr]*_lambda/np.pi
    print(the_arr_time)
    arr_times.append(the_arr_time)
# system2 = xy_chain.HeisenbergXY(t, N, 1, beta, corr='therm', lamb=_lambda)
# print(f"system: {time.time()-testing}")
# system.integrate()
# print(f"integrate {time.time()-testing}")
# system.i_dot_sq()
# print(f"inf_sq {time.time()-testing}")
# system.e_dot_test()
# print(f"edot_test {time.time()-testing}")
# system.single_state_rel_ent()
# print(f"rel_ent {time.time()-testing}")
# rel_ent = system.quantities['rel_ent']
# system2.single_state_rel_ent()
# print(f"rel_ent {time.time()-testing}")
# rel_ent2 = system2.quantities['rel_ent']
# edot = (system.quantities['e_dot_test'])*np.pi/3
# edot_2 = np.diff(system.integrated[-1])/dt
# idot = system.quantities['i_dot_sq']

fig, ax = plt.subplots(figsize=(8/2.54, 9/5.08))
ax.plot(poss_alpha, arr_times, color = cmap(.5))
fig.tight_layout()
fig.subplots_adjust(left=.2, bottom=.24)
ax.set_xlabel(r'$\alpha/\alpha_\mathrm{max}$')
ax.set_ylabel(r'Reduced Arrival Time')
plt.savefig('arrival_times.pdf', backend='pgf')

"""
ax.plot(t*_lambda, rel_ent, color = cmap(1/3), label = r'$\alpha=0$')
ax.plot(t*_lambda, rel_ent2, color = cmap(2/3), label = r'$\alpha=\alpha_\mathrm{max}$')
ax.set_yscale('log')
fig.tight_layout()
fig.subplots_adjust(left=.2, bottom=.24)
ax.set_xlabel(r'Time $\lambda t$')
ax.set_ylabel(r'Relative Entropy')
ax.legend(loc='best', frameon=False)
plt.savefig('state_transfer_corr_vs_uncorr.pdf', backend='pgf')

for i in range(N):
    _label = r'$\expval{\sigma^z_{'+str(i+1)+'}}$'
    ax.plot(t*_lambda, system.integrated[i], color=cmap((i+1)/(N+2)),
            label=_label)

_label = r'$\expval{\sigma^z_{'+str(N+1)+'}}$'
ax.plot(t*_lambda, system.integrated[-1], color=cmap(6/7), label = _label)
# ax.set_ylim(-1,.3)
# ax.legend(loc='best', ncol=2, frameon=False)
fig.tight_layout()
fig.subplots_adjust(left=.2,bottom=.24)
ax.set_xlabel(r'Time $\lambda t$')
ax.set_ylabel(r'Magnetization')
# ax.legend(loc='best')
plt.savefig('magnetization_no_corr.pdf', backend='pgf')

cmap = cm['viridis']

fig, ax = plt.subplots(figsize=(8/2.54,9/5.08))
_label = r'$\frac{\pi}{3}\dot{E}$'
ax.plot(t*_lambda, edot, color = cmap(1/3), label = _label)
_label = r'$\dot{I}^2$'
ax.plot(t[:-1]*_lambda, idot, color = cmap(2/3), label = _label)
ax.legend(loc='best', frameon=False)
fig.tight_layout()
fig.subplots_adjust(left=.2,bottom=.24)
ax.set_xlabel(r'Time $\lambda t$')
ax.set_ylabel(r'Information- and Energy flow', fontsize=7)
plt.savefig("pendry_no_corr.pdf", backend='pgf')
fig, ((ax0, ax1), (ax2,ax3)) = plt.subplots(2,2,figsize=(24,13.5))

for i in range(N):
    print(i)
    doft = system.n_qubit_discord_of_t([i,i+1])
    inds = np.argmax(doft)
    crit_rho[f'{i},{i+1}'] = system.time_evo[inds]
    print(crit_rho[f'{i},{i+1}'], crit_rho[f'{i},{i+1}'].ptrace([i,i+1]))
    print()
    print(crit_rho[f'{i},{i+1}'].eigenenergies(), crit_rho[f'{i},{i+1}'].ptrace([i,i+1]).eigenenergies())
    print()
    ax0.plot(t, doft, color=cmap((i+1)/(N+1)), label=r'$D(\rho_{'+f'{i+1},{i+2}'+r'})$')
    ax2.plot(t, system.integrated[i], color=cmap((i+1)/(N+2)))

ax2.plot(t, system.integrated[-1], color=cmap(6/7))

ax0.legend(loc='best')

ddot = np.diff(doft)/dt

ax1.plot(t[1:], edot*np.log(2)**2, ls=':', color = 'grey', label=r'no ln')
ax1.plot(t[1:], edot*np.log(2)**2+ddot**2, color = cmap(1/3), label=r'no ln + $(\partial_t D)^2$')
ax1.plot(t[1:], edot, color='grey', ls='-.', label=r'ln')
ax1.plot(t[1:], idot, color = cmap(2/3))
ax1.axhline(0, color='grey', ls='--')
ax1.legend(loc='best')

ax3.plot(t[1:], edot*np.log(2)**2-idot+ddot**2, color = cmap(1/2))
ax3.axhline(0, color='grey', ls='--')
fig.tight_layout()
plt.savefig("take_peek08.pdf")
plt.show()

#with open("data.pickle", "wb") as f:
#    pickle.dump(crit_rho, f)
# for ind, reduced_alpha in enumerate(poss_alpha):
#     time_list = [time.time(), time.time()]
#     system = xy_chain.HeisenbergXY(t, N, reduced_alpha, beta)
#     system.arrival_time()
#     arr_times.append(system.__dict__['poss_t'][0])
#     discord.append(system.discord)
#     time_list[1] = time.time()
#     print(ind, np.abs(np.diff(time_list)))
# xy_chain.write_alpha_checkpoint(
# [arr_times, discord],
# "pickled_data/N="+str(N+1)+"/arr_time_for_reduced_alpha="
# +str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)
# +"_"+str(poss_alpha.index(1)))
# xy_chain.write_alpha_checkpoint(
# [arr_times, discord],
# "pickled_data/N="+str(N+1)+"/arr_time_for_reduced_alpha"
# )
"""


