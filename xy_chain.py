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
# so far: only system with correlations at position 1&2, possibly subject to change
# 
# playing around with alpha and beta can give some arbitrary uncorrelated initial state at position 1


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


def init_system_arb_corr(first_state, nth_state, corr, alpha, N):
    qubit_states = [Qobj([[0,0],[0,1]]) for i in range(N+1)]
    qubit_states[corr[0]] = first_state
    qubit_states[corr[1]] = nth_state
    
    rho_system = qubit_states[0]    # calculate product state of system

    for state in qubit_states[1:]:
        rho_system = tensor(rho_system, state)

    # calculate longrange correlation term.
    # this part of the code is proper garbage, but mathematically correct
    # i need to think of something more efficient here
    identities = [qeye(2) for i in range(corr[1]-corr[0])]  
    sigma_p = [sigmap()]
    sigma_m = [sigmam()]
    sigma_p_0 = sigma_p.copy()
    sigma_p_0.extend(identities)
    sigma_m_0 = sigma_m.copy()
    sigma_m_0.extend(identities)
    sigma_p_1 = sigma_p.copy()
    sigma_p_1.extend(identities)
    sigma_p_1.reverse()
    sigma_m_1 = sigma_m.copy()
    sigma_m_1.extend(identities)
    sigma_m_1.reverse()

    chi = alpha * (tensor(sigma_p_0) * tensor(sigma_m_1) - tensor(sigma_m_0) * tensor(sigma_p_1))

    chi_helper = [qeye(2) for i in range(N+1)]
    chi_helper[corr[0]] = chi
    del chi_helper[corr[0]+1:corr[1]+1]
    total_chi = tensor(chi_helper)
    rho_system = rho_system + total_chi
    return rho_system


def thermal_state(beta):
    rho = (-beta * sigmaz()).expm()
    Z = rho.tr()
    return rho/Z


class System: 
    def __init__(self, t, N, H_given, alpha_reduced, beta: list, corr: list) -> None:
        self.t = t                      # array of timesteps
        self.dt = self.t[1]             # length of timestep
        self.N = N                      # Number of qubits
        self.H = H_given                # Hamiltonian
        self.beta = beta                # list of temperatures of the correlated qubits
        
        self.first_state = thermal_state(self.beta[0])  # first state in the correlation pair
        self.second_state = thermal_state(self.beta[1]) # second state in the correlation pair

        self.alpha = alpha_reduced * 1.j / (4*np.cosh(self.beta[0])*np.cosh(self.beta[1]))
        self.prod_state = tensor(self.first_state, self.second_state) # product state
        if type(corr) == None:
            #self.prod_state = tensor(self.first_state, self.second_state) # product state
            self.chi = self.alpha * tensor(Qobj([[0,1],[0,0]]), Qobj([[0,0],[1,0]])) - self.alpha * tensor(Qobj([[0,0],[1,0]]), Qobj([[0,1],[0,0]]))
            self.corr_state = self.prod_state + self.chi
            self.rho = init_system(self.alpha, self.first_state, self.second_state, self.N)
        else:    
            self.corr = corr                # list of qubits to correlate. if None, defaults to [0,1]
            self.corr_state = self.prod_state + self.alpha * tensor(Qobj([[0,1],[0,0]]), Qobj([[0,0],[1,0]])) - self.alpha * tensor(Qobj([[0,0],[1,0]]), Qobj([[0,1],[0,0]]))
            self.rho  = init_system_arb_corr(first_state=self.first_state, nth_state=self.second_state, corr=self.corr, alpha = self.alpha, N=self.N)

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

    
    def e_dot(self): # differentiating self.integrated does the same thing quicker (so far)
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