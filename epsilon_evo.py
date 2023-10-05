import pickle
import numpy as np
import time
from qutip import *
import matplotlib.pyplot as plt

import xy_chain


def epsilon_t(rho0, ham, eps):
    rho_small_t = rho0 + 1.j*eps*commutator(ham,rho0)
    return rho_small_t

cmap = plt.get_cmap('viridis')

x_state = Qobj([
    [.4, 0, 0, 1.j*.04],
    [0, .3, 1.j*.06, 0],
    [0, -1.j*.06, .2, 0],
    [-1.j*.04, 0, 0, .1]
])
#x_state = tensor(Qobj([[0, -.4], [-.6, 0]]), sigmay())

t_stop = np.pi
dt = 1e-3
t = np.arange(0, t_stop, dt)
beta = [1e-10, 1]
N = 1

system = xy_chain.HeisenbergXY(t, N, 1, beta, corr='therm')

rho_0 = system.init_system()

small_num = np.logspace(-10,-1,10, endpoint=True)
op_label = [["i", "x", "y", "z"]] 

sp1 = tensor(sigmap(), qeye(2))
sm1 = tensor(sigmam(), qeye(2))

sp2 = tensor(qeye(2), sigmap())
sm2 = tensor(qeye(2), sigmam())

x_state = .25*(tensor(qeye(2),qeye(2))
              +1.j*(sp1*sm2 + sp1*sp2 - sm1*sm2 - sm1*sp2))
print(x_state)
print(xy_chain.quantum_discord(x_state))
print(x_state.eigenenergies())
chi = system.chi
#rho_00 = tensor(system.uncorr,Qobj([[0, 0], [0, 1]]))
hkk1 = sp1*sm2+sm1*sp2-tensor(sigmaz(),sigmaz())
# Qobj([
#     [0, 0, 0, 0],
#     [0, 0, 1, 0],
#     [0, 1, 0, 0],
#     [0, 0, 0, 0]
# ])#
print(hkk1)
print(1.j*commutator(hkk1, x_state), 1.j*commutator(hkk1,chi))
showme = [rho_0, chi,
    #.01j*commutator(system.H, rho_00),
    .01j*commutator(system.H, chi),
    epsilon_t(rho_0, system.H, .01)]
print(system.H, chi, showme[2])
#print(.01j*commutator(system.H, chi))
#print(chi)
#print(rho_0+.01j*commutator(system.H, rho_0))
#system.integrate()

#fig, ax2 = plt.subplots(figsize=(12,5))
#for i in range(N+1):
#    ax2.plot(t, system.integrated[i], color=cmap((i+1)/(N+2)))
#plt.show()
#fig0, ax0 = matrix_histogram_complex(showme[0])
#fig1, ax1 = matrix_histogram_complex(showme[1])
#fig2, ax2 = matrix_histogram_complex(showme[2])
#fig3, ax3 = matrix_histogram_complex(showme[3])
#fig4, ax4 = matrix_histogram_complex(showme[4])
#for mat in showme:
#    print(mat)
##    print(epsilon_t(rho_0,system.H,e))
#    #fig, ax = 
#plt.show()
