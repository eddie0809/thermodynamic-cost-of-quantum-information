from unittest import result
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import linalg as LA
import prep_state
from time_evo import *
import my_plots
from generate_states import *
from qutip import *
import partialtrace

# for using tex formatting and font in plots

plt.rcParams.update({"text.usetex": True,}) 
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}\inputencoding{utf8}\usepackage{amsmath}\usepackage{amssymb}\usepackage{dsfont}\usepackage{mathtools}\usepackage{physics}'
mpl.rcParams['font.family'] = ['serif']

cmap = mpl.cm.get_cmap('viridis')

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
	for qubit, rho in enumerate(qubit_states):
		rho_system = tensor(rho_system, rho)
	return rho_system


t_stop = np.pi
dt=1e-3
t = np.arange(0,t_stop,dt)

beta = 1/300
rho = (-beta * sigmaz()).expm()
Z = rho.tr()
rho = rho/Z

N = 5
pos=2
alpha = 1.j * .25 * 1/(np.cosh(beta))**2

Hint = generate_hamiltonian(N)

first_state = Qobj([[1,0],[0,0]])
prod_state = tensor(rho, rho)
rho_system_test = system_arbitrary_pos_of_corr(alpha, first_state, prod_state, pos, N)

time_evolved = qutip.mesolve(Hint, rho_system_test, t, [], []).states
print(time_evolved[-1].ptrace(N))
fidelity_of_t = [fidelity(Qobj([[1,0],[0,0]]), timestep_evolved.ptrace(N)) for ind, timestep_evolved in enumerate(time_evolved)]


fig, (ax, axs1) = plt.subplots(2, 1, sharex=True, figsize=(12/2.54,9/2.54))
sz = integrate(N, rho_system_test, t, "me", "perf")

for n in range(N+1):
    ax.plot(t, sz[n], label=r"$\expval{\sigma^z_{%d}}$"%(n+1), lw=1.5, color=cmap((n+1)/(N+2)))
ax.set_xlabel(r'Time')
ax.set_ylabel(r'$\expval{\sigma^z_i}$')
axs1.plot(t, fidelity_of_t, label=r"$\alpha = \alpha_\mathrm{max}$", color=cmap(1/3))
fig.legend(loc=8, ncol=4)
fig.tight_layout()
fig.subplots_adjust(bottom=.26)
plt.show()

#print(rho_system_test, rho_system_test.tr())


"""

maxEvo = qutip.mesolve(Hint, rho_system, t, [], []).states



maxFid = [fidelity(first_state, maxEvo[dt].ptrace(N)) for dt in range(len(t))]


fig, axs1 = plt.subplots(4,1,sharex=True, sharey=True, figsize=(12,6))
axs1.plot(t, maxFid, label="alpha = alpha max")


fig.legend(loc=7)
fig.tight_layout()
fig.subplots_adjust(right=0.8)

plt.savefig("testingdingsplot.pdf")
plt.show()


# kullback leibler divergence

qEntmax = [entropy_vn(maxEvo[dt].ptrace(N), base=np.e) for dt in range(len(t))]
print('#')
infoFlowmax = np.square(np.diff(qEntmax)/t[1])
print('##')
maxHeatlim = [-1.j *(np.pi/3 *1 / (np.log(2)**2)) * (commutator(Hint, maxEvo[dt]).ptrace(N) * sigmaz()).tr() for dt in range(len(t))]
print('###')
fig, ax1 = plt.subplots(4,1,sharex=True, sharey=True, figsize=(12,6))
print('####')

ax1.plot(t, maxHeatlim, label="alpha = alpha max")
ax1.plot(t[:-1], infoFlowmax, label="alpha = alpha max")
print('-')

#ax1.axhline(0,color='grey')
#ax2.axhline(0,color='grey')
#ax3.axhline(0,color='grey')


fig.legend(loc=7)
fig.tight_layout()
fig.subplots_adjust(right=0.8)

plt.savefig("testplot.pdf")
plt.show()

"""
"""
def myqEnt(N):
	qEnt = [entropy_vn(timeEvo(dt, init_system, Hint).ptrace(N), base=np.e) for dt in t]
	myInfoflow = [(qEnt[i+1] - qEnt[i])/t[1] for i in range(0,len(qEnt)-1)]
	myHeatlim = [-1.j *(np.pi/3 *1 / (np.log(2)**2)) * (commutator(Hint, result[dt]).ptrace(N) * sigmaz()).tr() for dt in range(len(t))]
	#myHeatlim = [(np.pi/3 *1 / (np.log(2)**2)) * expect(timeEvo(dt, init_system, Hint).ptrace(N), sigmaz()) for dt in t]
	#myHeatlim = deriv(t[1], myHeatlim)
	#myHeat = [-1.j * (commutator(timeEvo(dt, init_system, Hint), Hint).ptrace(N) * -sigmaz()).tr() for dt in t]
	return qEnt, myInfoflow, myHeatlim#, myHeat
myKLdiv = [entropy_relative(first_state.unit(), timeEvo(dt, init_system, Hint).ptrace(N)) for dt in t]
#print(myKLdiv[-1])

helper = myqEnt(N)[1] 
myInfoflow2 = [helper[n]**2 for n in range(0,199)]

myTest = [-1.j * (commutator(timeEvo(dt, init_system, Hint), Hint).unit().ptrace(N) * sigmay()).tr() for dt in t]
#print(myTest)
myHeat = [-1.j * (commutator(timeEvo(dt, init_system, Hint), Hint).unit().ptrace(N) * sigmaz()).tr() for dt in t]
#print(np.mean(myHeat))
#print([commutator(result[dt], Hint).ptrace(N) for dt in range(len(t))])
fig, ax = plt.sub(figsize=(32/3,6))
#ax.plot(t, myKLdiv, label="$\mathcal{D}_\\text{KL}(\\rho^1(0)||\\rho^N(t))$")
#ax.set_title("Kullback-Leibler Divergence (relative entropy) \\\\ from $\\rho^1(0)$ to $\\rho^N(t)$ with $\\rho^1(0) = e^{-\\beta\sigma^z}/\Tr[e^{-\\beta\sigma^z}]$")
# 
# #\\frac{e^{-\\beta\sigma^z}}{\Tr[e^{-\\beta\sigma^z}]}

#plt.plot(t[:-1], myInfoflow2, label="$\dot{\mathcal{I}}_N^2$")
#plt.plot(t, myqEnt(N)[2], label="$\dot{E}_N\pi/(3\ln^2(2))$")
#ax.plot(t, myHeat, label="$\dot{\mathcal{Q}}_0$")
ax.set_xlabel("time")
#ax.set_ylabel("$\dot{\mathcal{Q}}_0$")
#ax.set_title("Heat flow at Qubit 0, $\dot{\mathcal{Q}}_0$. The grey bar indicates its limit predicted by Pendry, 1983.")
#plt.show()
#ax.plot(t[:-1], deriv(t[2]-t[1], myKLdiv), label="$\dot{\mathcal{D}}_\\text{KL}$")
#plt.title("Information flow from $\\rho^1(0)$ to $\\rho^N(t)$")
#plt.xlabel("time")
#plt.ylabel("$\dot{\mathcal{D}}_\\text{KL}$")
#plt.show()
#ax.plot(t,qEnt, label="$\mathcal{S}_{vn}$")
#ax.plot(t, np.real(myHeat), label="$\dot{\mathcal{Q}}$")
#ax.plot(t[:-1], dumbHeat)
ax.set_title("$\dot{\mathcal{I}}^2$ and its upper bound, $\dot{E}\pi/(3\ln^2(2))$ [Pendry, 1983], with $\\rho^0(0) = e^{-\\beta\sigma^z}/\Tr[e^{-\\beta\sigma^z}]$")
#plt.xlabel("time")
#plt.ylabel("$\dot{Q}_\\text{KL}$")
ax.legend()
#ax.axhline(np.pi/3 /(beta**2),color='grey')
ax.axhline(0,color='grey')
#plt.savefig("lost.png")
#plt.show()
plt.close('all')
"""