#!/usr/bin/env python3

import numpy as np
from numpy import exp, cosh, sqrt
from qutip import *
from q_discord import quantum_discord
import matplotlib as mpl
import matplotlib.pyplot as plt

# for using tex formatting and font in plots

plt.rcParams.update({"text.usetex": True,}) 
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}\inputencoding{utf8}\usepackage{amsmath}\usepackage{amssymb}\usepackage{dsfont}\usepackage{mathtools}\usepackage{physics}'
mpl.rcParams['font.family'] = ['serif']

cmap = mpl.cm.get_cmap('viridis')

steps = 100     # number of timesteps
D = np.zeros(steps)
a = np.linspace(0,1,steps)

b1 = 0
b2 = 0

p1 = exp(-b1)/(2*cosh(b1)) #np.array([[exp(-b1), 0], [0, exp(b1)]])
p2 = exp(-b2)/(2*cosh(b2)) #np.array([[exp(-b2), 0], [0, exp(b2)]]) 

A = Qobj([[p1, 0],[0, 1-p1]])
B = Qobj([[p2, 0],[0, 1-p2]])

alpha = sqrt((p1+p2-2*p1*p2)**2-(p1-p2)**2)/2

chi = [1j*a[j]*alpha*(tensor(basis(2,0),basis(2,1))*tensor(basis(2,1),basis(2,0)).dag() - tensor(basis(2,1),basis(2,0))*tensor(basis(2,0),basis(2,1)).dag()) for j in range(steps)]


for n in range(steps):
    rho = tensor(A,B) + chi[n]
    D[n] = quantum_discord(rho)

fig = plt.figure(figsize=(12/2.54,5/2.54))
plt.plot(a,D, color = cmap(1/2), label=r'Geometric Discord')
plt.xlabel(r'$\alpha/\alpha_\text{max}$')
plt.ylabel(r'$D(\alpha)$')
fig.tight_layout()
plt.savefig("geom_discord.pdf")
plt.show()
