import numpy as np
from generate_states import generate_hamiltonian

import xy_chain


t_stop = np.pi
dt = 1e-3
t = np.arange(0, t_stop, dt)
beta = [1/300, 1/300]
N = 5
Hint = generate_hamiltonian(N)

poss_alpha = [i/100 for i in range(101)]

for ind, reduced_alpha in enumerate(poss_alpha):
	system = xy_chain.System(t, N, Hint, reduced_alpha, beta)
	system.calculate_quantities()
	xy_chain.write_alpha_checkpoint(system.quantities, "pickled_data/N="+str(N+1)+"/quantities_with_reduced_alpha="+str(ind)+"_"+str(poss_alpha.index(1)))
	xy_chain.write_alpha_checkpoint(system.integrated, "pickled_data/N="+str(N+1)+"/integration_with_reduced_alpha="+str(ind)+"_"+str(poss_alpha.index(1)))
