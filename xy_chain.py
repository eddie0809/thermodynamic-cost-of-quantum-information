#!/usr/bin/env python3

import os
import pickle
import numpy as np

from qutip import *
from scipy.signal import argrelmin

from q_discord import quantum_discord


class HeisenbergXY:
    """There used to be two more functions.

    One returned a system with initial density matrix, where
    correlations could be at an arbitrary position.
    It can be found in a separate file named 'dumping_ground'.

    The other one was a function which computed "long-range correlations",
    which turned out to not work due to the resulting density matrix no longer
    being positive. check previous versions for it if you must.
    """
    def __init__(self, t, N, alpha_reduced,
                 beta: list, corr: str, x_state: Qobj = None, lamb = None) -> None:
        self.t = t             # array of timesteps
        self.dt = self.t[1]    # length of timestep
        self.N = N             # Number of qubits
        self.beta = beta       # list of temperatures of the correlated qubits
        self.alpha = 1. * alpha_reduced * self.get_max_alpha()

        if corr == "therm" and alpha_reduced > 0:
            self.first_state = (-beta[0]*sigmaz()).expm()/(2*np.cosh(beta[0]))
            second_state = (-beta[1]*sigmaz()).expm()/(2*np.cosh(beta[1]))
            self.uncorr = tensor(self.first_state, second_state)

            self.chi = self.alpha * tensor(
                 sigmap(), sigmam()) - self.alpha * tensor(sigmam(), sigmap())
            self.corr_state = self.uncorr+self.chi

            self.rho = self.init_system()

            self.discord = quantum_discord(self.corr_state)
        elif corr == "X":
            self.corr_state = x_state
            self.rho = self.init_system()
            # self.discord = quantum_discord(self.corr_state)
        elif alpha_reduced == 0 and corr != "X":
            self.first_state = (-beta[0]*sigmaz()).expm()/(2*np.cosh(beta[0]))
            second_state = (-beta[1]*sigmaz()).expm()/(2*np.cosh(beta[1]))
            self.corr_state = tensor(self.first_state, second_state)
            self.rho = self.init_system()

        self.calc_composite_ops()
        self.compute_xy_hamiltonian(lamb)

        self.time_evo = mesolve(self.H, self.rho, self.t, [], []).states
        self.quantities = {}

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

    def compute_xy_hamiltonian(self,  lamb=1., coupling=1/2):
        hamiltonian = 0
        for i in range(self.N):
            hamiltonian += (coupling*lamb) / 2 * np.sqrt((i+1)*(self.N-i)) * (
                 self.sx_list[i] * self.sx_list[i+1]
                 + self.sy_list[i] * self.sy_list[i+1])
            hamiltonian -= self.sz_list[i]
        self.H = hamiltonian - self.sz_list[self.N]

    def get_max_alpha(self):
        return 1.j / (4*np.cosh(self.beta[0])*np.cosh(self.beta[1]))

    def init_system(self):
        rho_system = self.corr_state
        for _ in range(1, self.N):
            rho_system = tensor(rho_system, Qobj([[0, 0], [0, 1]]))
        return rho_system

    def n_qubit_discord_of_t(self, qubits: list):
        return [quantum_discord(rho_t.ptrace(qubits)) for rho_t in
                self.time_evo]

    def single_state_fidelity(self):
        self.quantities["fidelity"] = np.array(
             [fidelity(self.first_state, time_evolved.ptrace(self.N))
              for time_evolved in self.time_evo])
        self.single_fidelity = self.quantities["fidelity"]
        self.quantities["fid_max"] = np.amax(self.single_fidelity)

    def single_state_rel_ent(self):
        self.quantities["rel_ent"] = np.array(
             [entropy_relative(self.first_state, time_evolved.ptrace(self.N))
              for time_evolved in self.time_evo])
        self.single_rel_ent = self.quantities["rel_ent"]

    def arrival_time(self):
        # self.single_state_fidelity()
        self.single_state_rel_ent()
        self.poss_t = argrelmin(self.single_rel_ent)
        self.quantities["local minima in rel ent"] = self.poss_t

    def i_dot_sq(self):
        vn_ent = [entropy_vn(time_evolved.ptrace(self.N))
                  for time_evolved in self.time_evo]
        self.quantities["i_dot_sq"] = np.square(np.diff(vn_ent)/self.dt)
        self.quantities["i_dot_sq_max"] = np.max(self.quantities["i_dot_sq"])

    def e_dot(self):
        # differentiating self.integrated does the same thing quicker (so far)
        self.quantities["e_dot"] = [
              -1.j*(commutator(self.H, evolved).ptrace(self.N)*sigmaz()).tr()
              for evolved in self.time_evo]

    def e_dot_test(self, k=None):
        # explainer: in Lie algebras we have tr(x * [y,z]) = tr([x,y] * z)
        # further, the von Neumann equation for the time derivative of rho
        # reads d/dt rho = -i [H, rho], i.e.
        # d/dt <x> = d/dt tr(x * rho) = -i * tr(x * [H, rho])
        # = - Im{tr([x,H] * rho)}
        # by somehow losing the - in the process, we can save on computational
        # cost by calculating the commutator once and skip painstakingly
        # computing the partial trace _and_ the commutator at each timestep
        # and instead compute one commutator and then only multiply matrices
        # and trace it.
        if k is None: # this makes sure that edot of N is computed by default
            k=self.N
        this_chi = commutator(self.sz_list[k], self.H)
        e_dot_of_t_test = np.imag([(this_chi*evolved).tr() for evolved in
                           self.time_evo])
        self.quantities["e_dot_test"] = e_dot_of_t_test

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
