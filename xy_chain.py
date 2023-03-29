#!/usr/bin/env python3

import os
import pickle
import numpy as np

from scipy.signal import argrelmin
from qutip import *

from q_discord import quantum_discord

# goal: system(10, 3, H, [a,b]) initializes a spinchain with given length (10), position of left qubit in corr state (3), an interaction H and two temperatures
# for the thermal states to have, all at time t=0.
#
# system(N, pos_corr, H_given, beta).time_evo(t) should return the time evolution for some time t (numpy array)
#
# so far: only system with correlations at position 1&2, possibly subject to change
# 
# playing around with alpha and beta can give some arbitrary uncorrelated initial state at position 1

"""
There used to be two functions here.

One returned a density matrix, where
correlations could be at an arbitrary position. It can be found in a separate
file named 'dumping_ground'.

The other one was a function which computed "long-range correlations",
which turned out to not work due to the resulting density matrix no longer
being positive. check previous versions for it if you must.
"""

class HeisenbergXY: 
    def __init__(self, t, N, alpha_reduced, beta: list) -> None: # , corr: list
        self.t = t                      # array of timesteps
        self.dt = self.t[1]             # length of timestep
        self.N = N                      # Number of qubits
        self.beta = beta                # list of temperatures of the correlated qubits
        self.alpha = alpha_reduced * self.get_max_alpha()

        self.first_state  = (-beta[0]*sigmaz()).expm()/(2*np.cosh(beta[0]))
        second_state = (-beta[1]*sigmaz()).expm()/(2*np.cosh(beta[1]))
        
        self.chi = self.alpha * tensor(sigmap(), sigmam()) - self.alpha * tensor(sigmam(), sigmap())
        self.corr_state = tensor(self.first_state, second_state) + self.chi
        
        self.rho = self.init_system()

        self.discord = quantum_discord(self.corr_state)
        
        self.calc_composite_ops()
        self.compute_xy_hamiltonian()

        self.time_evo = mesolve(self.H, self.rho, self.t, [], []).states
        self.quantities = {}
        self.meta = {"t":self.t, "dt":self.dt,"N":self.N, "beta":beta, "H":self.H, "discord":self.discord}


    def calc_composite_ops(self):
        self.sx_list = []
        self.sy_list = []
        self.sz_list = []
        self.sp_list = []
        self.sm_list = []

        for n in range(self.N+1):
            op_list = [qeye(2) for _ in range(self.N+1)]

            op_list[n] = sigmax()
            self.sx_list.append(tensor(op_list))

            op_list[n] = sigmay()
            self.sy_list.append(tensor(op_list))

            op_list[n] = sigmaz()
            self.sz_list.append(tensor(op_list))

            op_list[n] = sigmap()
            self.sp_list.append(tensor(op_list))

            op_list[n] = sigmam()
            self.sm_list.append(tensor(op_list))


    def compute_xy_hamiltonian(self,  l=1., coupling=1/2):
        hamiltonian = 0
        for i in range(self.N):
            hamiltonian += (coupling*l) / 2 * np.sqrt((i+1)*(self.N-i)) * (self.sx_list[i] * self.sx_list[i+1] + self.sy_list[i] * self.sy_list[i+1])
            hamiltonian -= self.sz_list[i]
        self.H = hamiltonian - self.sz_list[self.N]


    def get_max_alpha(self):
        return 1.j / (4*np.cosh(self.beta[0])*np.cosh(self.beta[1]))


    def init_system(self):
        rho_system = self.corr_state
        for _ in range(1,self.N):
            rho_system = tensor(rho_system, Qobj([[0,0],[0,1]]))
        return rho_system


    def single_state_fidelity(self):
        self.quantities["fidelity"] = np.array([fidelity(self.first_state, time_evolved.ptrace(self.N)) for time_evolved in self.time_evo])
        self.single_fidelity = self.quantities["fidelity"]
        self.quantities["fid_max"] = np.amax(self.single_fidelity)
    

    def single_state_rel_ent(self):
        self.quantities["rel_ent"] = np.array([entropy_relative(self.first_state, time_evolved.ptrace(self.N)) for time_evolved in self.time_evo])
        self.single_rel_ent = self.quantities["rel_ent"]


    def arrival_time(self):
        self.single_state_fidelity()
        self.single_state_rel_ent()
        poss_t = argrelmin(self.single_rel_ent)
        self.quantities["local minima in rel ent"] = poss_t
        

    def i_dot_sq(self):
        vn_ent = [entropy_vn(time_evolved.ptrace(self.N)) for time_evolved in self.time_evo]
        self.quantities["i_dot_sq"] = np.square(np.diff(vn_ent)/self.dt)
        self.quantities["i_dot_sq_max"] = np.max(self.quantities["i_dot_sq"])

    
    def e_dot(self): # differentiating self.integrated does the same thing quicker (so far)
        self.quantities["e_dot"] = [-1.j*(commutator(self.H, time_evolved).ptrace(self.N) * (sigmaz())).tr() for time_evolved in self.time_evo]


    def integrate(self, solver="me"):
        # evolve and calculate expectation values
        if solver == "me":
            result = mesolve(self.H, self.rho, self.t, [], self.sz_list)
        elif solver == "mc":
            ntraj = 250 
            result = mcsolve(self.H, self.rho, self.t, [], self.sz_list, ntraj)
        self.integrated = result.expect


def write_alpha_checkpoint(system_quantities, path, overwrite=False):
	if os.path.exists(path) and not overwrite:
		raise RuntimeError("Data for this alpha already exist")
	with open(path, 'wb') as fp:
		pickle.dump(system_quantities, fp)