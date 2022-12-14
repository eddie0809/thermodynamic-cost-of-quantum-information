#!/usr/bin/env python3

import numpy as np
from numpy import exp, cosh, sqrt
from qutip import *
from q_discord import quantum_discord
import matplotlib.pyplot as plt

steps = 100     # number of timesteps
D = np.zeros(steps)
a = np.linspace(0,1,steps)

b1 = 10
b2 = 10

p1 = exp(-b1)/(2*cosh(b1)) #np.array([[exp(-b1), 0], [0, exp(b1)]])
p2 = exp(-b2)/(2*cosh(b2)) #np.array([[exp(-b2), 0], [0, exp(b2)]]) 

A = Qobj([[p1, 0],[0, 1-p1]])
B = Qobj([[p2, 0],[0, 1-p2]])

alpha = sqrt((p1+p2-2*p1*p2)**2-(p1-p2)**2)/2

chi = [1j*a[j]*alpha*(tensor(basis(2,0),basis(2,1))*tensor(basis(2,1),basis(2,0)).dag() - tensor(basis(2,1),basis(2,0))*tensor(basis(2,0),basis(2,1)).dag()) for j in range(steps)]


for n in range(steps):
    rho = tensor(A,B) + chi[n]
    D[n] = quantum_discord(rho)

plt.plot(a,D)
plt.show()
