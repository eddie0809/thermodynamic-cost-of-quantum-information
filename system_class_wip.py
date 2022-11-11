import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import prepState
from numpy import linalg as LA
from qutip import *

# goal: system(10, 3, H, [a,b]) initializes a spinchain with given length (10), position of left qubit in corr state (3), an interaction H and two temperatures
# for the thermal states to have, all at time t=0.
#
# system(N, pos_corr, H_given, beta).time_evo(t) should return the time evolution for some time t (numpy array)
#
# 

class system: 
    def __init__(self, N, pos_corr, H_given, beta) -> None:
        self.N = N
        self.pos_corr = pos_corr        # position of left qubit in correlation state
        self.H = H_given 
        self.beta = beta                # list of temperatures

        thermal_states = []
        alpha = 1.j / (4*np.cosh(self.beta[0])*np.cosh(self.beta[1]))
        chi = alpha * tensor(Qobj([[0,1],[0,0]]), Qobj([[0,0],[1,0]])) - alpha * tensor(Qobj([[0,0],[1,0]]), Qobj([[0,1],[0,0]]))
        corr_state = tensor(thermal_states[0], thermal_states[1]) + chi



    def time_evo(self, t):
        return qutip.mesolve(self.H, self, t, [], []).states