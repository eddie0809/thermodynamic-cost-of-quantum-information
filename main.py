import numpy as np
from qutip import *
import xy_chain

t_stop = 2*np.pi
dt = 1e-3
t = np.arange(0, t_stop, dt)
beta = [0, 0]
N = 5

poss_alpha = [i/12 for i in range(0,13)]

for ind, reduced_alpha in enumerate(poss_alpha):
	system = xy_chain.HeisenbergXY(t, N, reduced_alpha, beta)
	two_qubit_fidelity = [fidelity(system.corr_state, time_evolved_system.ptrace([N,N-1])) for time_evolved_system in system.time_evo]
	xy_chain.write_alpha_checkpoint(system.meta, "pickled_data/N="+str(N+1)+"/metadata_for_reduced_alpha="+str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)+"_"+str(poss_alpha.index(1)))
	xy_chain.write_alpha_checkpoint(system.integrated, "pickled_data/N="+str(N+1)+"/integration_with_reduced_alpha="+str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)+"_"+str(poss_alpha.index(1)))
	xy_chain.write_alpha_checkpoint(system.time_evo, "pickled_data/N="+str(N+1)+"/system_integrated_reduced_alpha="+str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)+"_"+str(poss_alpha.index(1)))
	xy_chain.write_alpha_checkpoint(two_qubit_fidelity, "pickled_data/N="+str(N+1)+"/new_two_qubit_fidelity_with_reduced_alpha="+str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)+"_"+str(poss_alpha.index(1)))
