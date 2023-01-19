import numpy as np
from generate_states import generate_hamiltonian
from qutip import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import xy_chain

# for using tex formatting and font in plots

plt.rcParams.update({"text.usetex": True,}) 
plt.rcParams['text.latex.preamble'] = r'\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}\inputencoding{utf8}\usepackage{amsmath}\usepackage{amssymb}\usepackage{dsfont}\usepackage{mathtools}\usepackage{physics}'
mpl.rcParams['font.family'] = ['serif']

cmap = mpl.cm.get_cmap('PuOr')

t_stop = np.pi
dt = 1e-3
t = np.arange(0, t_stop, dt)
beta = [0, .5]
N = 5
Hint = generate_hamiltonian(N)

poss_alpha = [0,.1,1]

for ind, reduced_alpha in enumerate(poss_alpha):
	system = xy_chain.System(t, N, Hint, reduced_alpha, beta)
	#print(system.time_evo[-1].ptrace([N-1,N]))
	print(system.rho.eigenenergies())
	plt.figure(figsize=(12,5))
	for ind, sz in enumerate(system.integrated):
		plt.plot(t, sz, label = str(ind), color = cmap((ind+1)/(N+2)))
	plt.legend()
	plt.show()
	#prod_state = tensor(system.second_state, system.first_state)
	#print(fidelity(prod_state-system.chi, system.time_evo[-1].ptrace([N-1,N])))
	#system.single_state_fidelity()
	#plt.figure(figsize=(14.5,8))
	#plt.plot(t, system.quantities['fidelity'])
	#plt.show()
	#two_qubit_fidelity = [fidelity(system.corr_state, time_evolved.ptrace([N,N-1])) for time_evolved in system.time_evo]
	#xy_chain.write_alpha_checkpoint(system.quantities, "pickled_data/N="+str(N+1)+"/quantities_with_reduced_alpha="+str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)+"_"+str(poss_alpha.index(1)))
	#xy_chain.write_alpha_checkpoint(system.integrated, "pickled_data/N="+str(N+1)+"/integration_with_reduced_alpha="+str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)+"_"+str(poss_alpha.index(1)))
	#xy_chain.write_alpha_checkpoint(system.time_evo, )
	#xy_chain.write_alpha_checkpoint(two_qubit_fidelity, "pickled_data/N="+str(N+1)+"/new_two_qubit_fidelity_with_reduced_alpha="+str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)+"_"+str(poss_alpha.index(1)))

