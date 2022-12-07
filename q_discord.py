#!/usr/bin/env python3

import numpy as np
from qutip import *

def quantum_discord(rho):
    x = np.zeros(3)*1j
    V = np.zeros((3,3))*1j
    sigma = np.array([sigmax(),sigmay(),sigmaz()])


    for i, sig in enumerate(sigma):
        x[i] = (rho*tensor(Qobj(sig),qeye(2))).tr() 
        for j, sigm in enumerate(sigma):
            V[i,j] = (rho*tensor(Qobj(sig),Qobj(sigm))).tr()

    x = Qobj(x)
    V = Qobj(V)
    Lambda = x*x.dag() + V*V.dag()

    k = max(Lambda.eigenenergies())

    D = 2*(Lambda.tr()-k)

    return D
