import numpy as np
from generate_states import generate_hamiltonian
from qutip import *
import matplotlib.pyplot as plt
import xy_chain


t_stop = np.pi
dt = 1e-3
t = np.arange(0, t_stop, dt)
beta = [0, 0]
N = 5
Hint = generate_hamiltonian(N)

poss_alpha = [0,.5,1]

for ind, reduced_alpha in enumerate(poss_alpha):
	system = xy_chain.System(t, N, Hint, reduced_alpha, beta, corr=[0,2])
	#print(system.time_evo[-1].ptrace([N-1,N]))
	#print(system.corr_state)
	for ind, sz in enumerate(system.integrated):
		plt.plot(t, sz, label = str(ind))
	plt.legend()
	plt.show()
	#prod_state = tensor(system.second_state, system.first_state)
	#print(fidelity(prod_state-system.chi, system.time_evo[-1].ptrace([N-1,N])))
	system.calculate_quantities()
	plt.plot(t, system.quantities['fidelity'])
	plt.show()
	#two_qubit_fidelity = [fidelity(system.corr_state, time_evolved.ptrace([N,N-1])) for time_evolved in system.time_evo]
	#xy_chain.write_alpha_checkpoint(system.quantities, "pickled_data/N="+str(N+1)+"/quantities_with_reduced_alpha="+str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)+"_"+str(poss_alpha.index(1)))
	#xy_chain.write_alpha_checkpoint(system.integrated, "pickled_data/N="+str(N+1)+"/integration_with_reduced_alpha="+str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)+"_"+str(poss_alpha.index(1)))
	#xy_chain.write_alpha_checkpoint(system.time_evo, )
	#xy_chain.write_alpha_checkpoint(two_qubit_fidelity, "pickled_data/N="+str(N+1)+"/new_two_qubit_fidelity_with_reduced_alpha="+str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)+"_"+str(poss_alpha.index(1)))

