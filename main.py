import numpy as np
from generate_states import generate_hamiltonian
from qutip import *
import xy_chain


t_stop = np.pi
dt = 1e-3
t = np.arange(0, t_stop, dt)
beta = [1/300, 1/300]
N = 5
Hint = generate_hamiltonian(N)

poss_alpha = [i/12 for i in range(13)]

for ind, reduced_alpha in enumerate(poss_alpha):
	system = xy_chain.System(t, N, Hint, reduced_alpha, beta)
	system.calculate_quantities()
	two_qubit_fidelity = [fidelity(system.corr_state, time_evolved.ptrace([N-1,N])) for time_evolved in system.time_evo]
	xy_chain.write_alpha_checkpoint(system.quantities, "pickled_data/N="+str(N+1)+"/quantities_with_reduced_alpha="+str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)+"_"+str(poss_alpha.index(1)))
	xy_chain.write_alpha_checkpoint(system.integrated, "pickled_data/N="+str(N+1)+"/integration_with_reduced_alpha="+str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)+"_"+str(poss_alpha.index(1)))
	xy_chain.write_alpha_checkpoint(two_qubit_fidelity, "pickled_data/N="+str(N+1)+"/two_qubit_fidelity_with_reduced_alpha="+str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)+"_"+str(poss_alpha.index(1)))