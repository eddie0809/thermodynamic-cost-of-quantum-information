#!/usr/bin/env python3

import os
import pickle

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import prep_state
from numpy import linalg as LA
from qutip import *
from q_discord import quantum_discord

# goal: system(10, 3, H, [a,b]) initializes a spinchain with given length (10), position of left qubit in corr state (3), an interaction H and two temperatures
# for the thermal states to have, all at time t=0.
#
# system(N, pos_corr, H_given, beta).time_evo(t) should return the time evolution for some time t (numpy array)
#
# test

def integrate(N, H, psi0, tlist, solver):

    si = qeye(2)
    sx = sigmax()
    sy = sigmay()
    sz = sigmaz()

    sx_list = []
    sy_list = []
    sz_list = []

    for n in range(N+1):
        op_list = [si for m in range(N+1)]

        op_list[n] = sx
        sx_list.append(tensor(op_list))

        op_list[n] = sy
        sy_list.append(tensor(op_list))

        op_list[n] = sz
        sz_list.append(tensor(op_list))

    # construct the hamiltonian
    #H = generate_states(N, interaction)

    # energy splitting terms

    # interaction terms
    

    # collapse operators
    c_op_list = []


    # evolve and calculate expectation values
    if solver == "me":
        result = mesolve(H, psi0, tlist, c_op_list, sz_list)
    elif solver == "mc":
        ntraj = 250 
        result = mcsolve(H, psi0, tlist, c_op_list, sz_list, ntraj)
    return result.expect


def init_system(alpha, first_state, second_state, N):
	astar = -alpha
	rho_system = tensor(first_state, second_state)+ alpha * tensor(Qobj([[0,1],[0,0]]), Qobj([[0,0],[1,0]])) + astar * tensor(Qobj([[0,0],[1,0]]), Qobj([[0,1],[0,0]]))
	for i in range(1,N):
		rho_system = tensor(rho_system, Qobj([[0,0],[0,1]]))
	return rho_system


def system_arbitrary_pos_of_corr(alpha, first_state, prod_state, pos, N):
	if pos == 1:
		return init_system(alpha, prod_state.ptrace(0),prod_state.ptrace(1), N)
	qubit_states = []
	i = 1
	astar = -alpha
	corr_state = prod_state + alpha * tensor(Qobj([[0,1],[0,0]]), Qobj([[0,0],[1,0]])) + astar * tensor(Qobj([[0,0],[1,0]]), Qobj([[0,1],[0,0]]))
	rho_system = first_state
	while i < pos-1:
		qubit_states.append(Qobj([[0,0],[0,1]]))
		i+=1
	qubit_states.append(corr_state)
	i+=2
	while i<=N:
		qubit_states.append(Qobj([[0,0],[0,1]]))
		i+=1
	for rho in qubit_states:
		rho_system = tensor(rho_system, rho)
	return rho_system


def thermal_state(beta):
    rho = (-beta * sigmaz()).expm()
    Z = rho.tr()
    return rho/Z


class Simulation: 
    def __init__(self, t, N, H_given, alpha_reduced, beta) -> None:
        self.t = t                      # array of timesteps
        self.N = N
        self.H = H_given                # Hamiltonian
        self.beta = beta                # list of temperatures
        self.dt = self.t[1]             # timestep
        self.alpha = alpha_reduced * 1.j / (4*np.cosh(self.beta[0])*np.cosh(self.beta[1]))
        
        self.first_state = thermal_state(self.beta[0])
        self.second_state = thermal_state(self.beta[1])

        prod_state = tensor(self.first_state, self.second_state)
        chi = self.alpha * tensor(Qobj([[0,1],[0,0]]), Qobj([[0,0],[1,0]])) - self.alpha * tensor(Qobj([[0,0],[1,0]]), Qobj([[0,1],[0,0]]))

        self.corr_state = prod_state + chi

        self.rho = init_system(self.alpha, self.first_state, self.second_state, self.N)

        self.discord = quantum_discord(self.corr_state)

        self.time_evo = qutip.mesolve(self.H, self.rho, self.t, [], []).states
        self.integrated = integrate(self.N, self.H, self.rho, self.t, "me")
        self.quantities = {"fidelity":[], "rel_ent":[], "i_dot_sq":[], "i_dot_sq_max":[], "e_dot":[], "time_where_kld_isclose_0":[], "fid_max":[], "discord":self.discord}

    
    def single_state_fidelity(self):
        self.quantities["fidelity"] = np.array([fidelity(self.first_state, time_evolved.ptrace(self.N)) for time_evolved in self.time_evo])
        self._fidelity = self.quantities["fidelity"]
    

    def single_state_rel_ent(self):
        self.quantities["rel_ent"] = np.array([entropy_relative(self.first_state, time_evolved.ptrace(self.N)) for time_evolved in self.time_evo])
        self._rel_ent = self.quantities["rel_ent"]


    def arrival_time(self):
        self.single_state_fidelity()
        self.single_state_rel_ent()
        poss_t = np.where(np.isclose(self._rel_ent, 0, atol=1e-8), self.t, np.nan)
        self.quantities["fid_max"] = np.amax(self._fidelity)
        self.quantities["time_where_kld_isclose_0"] = poss_t
        

    def i_dot_sq(self):
        vn_ent = [entropy_vn(time_evolved.ptrace(self.N)) for time_evolved in self.time_evo]
        self.quantities["i_dot_sq"] = np.square(np.diff(vn_ent)/self.dt)
    

    def i_dot_sq_max(self):
        self.i_dot_sq()
        self.quantities["i_dot_sq_max"] = np.max(self.quantities["i_dot_sq"])

    
    def e_dot(self):
        self.quantities["e_dot"] = [-1.j*(commutator(self.H, time_evolved).ptrace(self.N) * (sigmaz())).tr() for time_evolved in self.time_evo]


    def calculate_quantities(self):
        self.arrival_time()
        self.i_dot_sq_max()
        self.e_dot()


def write_alpha_checkpoint(system_quantities, path, overwrite=False):
	if os.path.exists(path) and not overwrite:
		raise RuntimeError("Data for this alpha already exist")
	with open(path, 'wb') as fp:
		pickle.dump(system_quantities, fp)