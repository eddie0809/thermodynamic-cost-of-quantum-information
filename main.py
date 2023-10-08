import pickle
import numpy as np
import time
from qutip import *
import matplotlib.pyplot as plt
from matplotlib import colormaps as cm

import xy_chain

cmap = cm['viridis']

beta = [1e-3, .01]
N = 7
t_stop = np.pi*np.sqrt(N)/2
dt = 1e-3
t = np.arange(0, t_stop, dt)
testing = time.time()
print(time.time()-testing)
x_state = Qobj([
    [.4, 0, 0, 1.j*.04],
    [0, .3, 1.j*.06, 0],
    [0, -1.j*.06, .2, 0],
    [-1.j*.04, 0, 0, .1]
])

poss_alpha = [i/100 for i in range(101)]

arr_times = []
discord = []
crit_rho = {}
system = xy_chain.HeisenbergXY(t, N, 1, beta, corr='therm')
print(f"system: {time.time()-testing}")
# system.integrate()
# print(f"integrate {time.time()-testing}")
system.i_dot_sq()
print(f"inf_sq {time.time()-testing}")
system.e_dot_test()
print(f"edot_test {time.time()-testing}")
system.single_state_rel_ent()
print(f"rel_ent {time.time()-testing}")
rel_ent = system.quantities['rel_ent']
edot = (system.quantities['e_dot_test'])*np.pi/3
# edot_2 = np.diff(system.integrated[-1])/dt
idot = system.quantities['i_dot_sq']

plt.plot(edot)
plt.plot(idot)
plt.plot(edot[:-1]-idot)
plt.show()
plt.plot(rel_ent)
plt.yscale('log')
# plt.xscale('log')
plt.show()
"""
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

