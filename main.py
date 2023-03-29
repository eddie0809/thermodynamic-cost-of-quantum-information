import numpy as np
import xy_chain

t_stop = 2*np.pi
dt = 1e-3
t = np.arange(0, t_stop, dt)
beta = [0, 0]
N = 5

poss_alpha = [i/12 for i in range(0,13)]

for ind, reduced_alpha in enumerate(poss_alpha):
	system = xy_chain.HeisenbergXY(t, N, reduced_alpha, beta)
	system.arrival_time()
	system.i_dot_sq()
	system.integrate()
	xy_chain.write_alpha_checkpoint(system.__dict__, "pickled_data/N="+str(N+1)+"/dict_for_reduced_alpha="+str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)+"_"+str(poss_alpha.index(1)))
	