import numpy as np
from qutip import *

import prep_state


def generate_hamiltonian(N):
    my_states = {}
    my_basis = {}
    for i in range(N+1):
        my_states["sigmaX", i] = prep_state.stateX(N,i)
        my_states["sigmaY", i] = prep_state.stateY(N,i)
        my_states["sigmaZ", i] = prep_state.stateZ(N,i)
        my_basis["zero",i] = prep_state.state0(N,i)
        my_basis["one",i] = prep_state.state1(N,i)

    l = 1
    Hint = 0
    for i in range(N): # interaction hamiltonian, XY-Heisenberg
        Hint = Hint + l / 2 * np.sqrt((i+1)*(N-i)) * (my_states["sigmaX", i] * my_states["sigmaX", i+1] + my_states["sigmaY", i] * my_states["sigmaY", i+1])
    J = 1/2
    Hint = J * Hint # change interaction strength 
    H0 = 0
    for i in range(0,N+1):
        H0 = H0 - my_states["sigmaZ", i]
    return H0 + Hint