import numpy as np
import xy_chain
import time

t_stop = np.pi
dt = 1e-3
t = np.arange(0, t_stop, dt)
beta = [0, 0]
N = 5

poss_alpha = [i/100 for i in range(101)]

arr_times = []
discord = []
for ind, reduced_alpha in enumerate(poss_alpha):
	time_list = [time.time(), time.time()]
	system = xy_chain.HeisenbergXY(t, N, reduced_alpha, beta)
	system.arrival_time()
	arr_times.append(system.__dict__['poss_t'][0])
	discord.append(system.discord)
	time_list[1] = time.time()
	print(ind, np.abs(np.diff(time_list)))
#xy_chain.write_alpha_checkpoint([arr_times, discord], "pickled_data/N="+str(N+1)+"/arr_time_for_reduced_alpha="+str(int((ind/100))%10)+str(int((ind/10))%10)+str(ind%10)+"_"+str(poss_alpha.index(1)))
xy_chain.write_alpha_checkpoint([arr_times, discord], "pickled_data/N="+str(N+1)+"/arr_time_for_reduced_alpha")	
