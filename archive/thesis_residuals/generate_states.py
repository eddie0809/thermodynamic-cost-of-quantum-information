#!/usr/bin/env python3

import numpy as np
from qutip import *


def calc_composite_ops(N):
    sx_list = []
    sy_list = []
    sz_list = []
    sp_list = []
    sm_list = []

    for n in range(N+1):
        op_list = [qeye(2) for _ in range(N+1)]

        op_list[n] = sigmax()
        sx_list.append(tensor(op_list))

        op_list[n] = sigmay()
        sy_list.append(tensor(op_list))

        op_list[n] = sigmaz()
        sz_list.append(tensor(op_list))

        op_list[n] = sigmap()
        sp_list.append(tensor(op_list))

        op_list[n] = sigmam()
        sm_list.append(tensor(op_list))
    return sx_list, sy_list, sz_list


def generate_xy_hamiltonian(N, l=1., coupling=1/2) -> Qobj:
    sx_list, sy_list, sz_list = calc_composite_ops(N)
    
    hamiltonian = 0

    for i in range(N):
        hamiltonian += l / 2 * np.sqrt((i+1)*(N-i)) * (sx_list[i] * sx_list[i+1] + sy_list[i] * sy_list[i+1] )
    
    hamiltonian *= coupling
    for i in range(0,N+1):
        hamiltonian -= sz_list[i]

    return hamiltonian