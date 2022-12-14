### this file is just that the main file isnt filled with 2 line functions.
import numpy as np
from qutip import *
from scipy import linalg as la

import prep_state


def generate_states(N, interaction):
    myStates = {}
    myBasis = {}
    for i in range(N+1):
        myStates["sigmaX", i] = prep_state.stateX(N,i)
        myStates["sigmaY", i] = prep_state.stateY(N,i)
        myStates["sigmaZ", i] = prep_state.stateZ(N,i)
        myStates["sigmaP", i] = .5*(myStates["sigmaX", i]+1.j*myStates["sigmaY", i])
        myStates["sigmaM", i] = .5*(myStates["sigmaX", i]-1.j*myStates["sigmaY", i])
        myBasis["zero",i] = prep_state.state0(N,i)
        myBasis["one",i] = prep_state.state1(N,i)

    #my_one = tensor(qeye(2),qeye(2),qeye(2),qeye(2),qeye(2))
    if interaction=="hom": #homogenous interaction strength. 
        Hint = 0
        for i in range(N): # interaction hamiltonian, XY-Heisenberg #l / 2 * np.sqrt((i+1)*(N-i))
            Hint = Hint + (myStates["sigmaX", i] * myStates["sigmaX", i+1] + myStates["sigmaY", i] * myStates["sigmaY", i+1])
        J = 1
        Hint = J * Hint # change interaction strength 

        H0 = 0
        for i in range(0,N+1):
                H0 = H0 + myStates["sigmaZ", i]
        Hges = Hint + H0
        return Hges
        
    elif interaction=="perf":
        Hint = 0
        l=1
        for i in range(N): # interaction hamiltonian, XY-Heisenberg #
            Hint = Hint + l / 2 * np.sqrt((i+1)*(N-i)) * (myStates["sigmaX", i] * myStates["sigmaX", i+1] + myStates["sigmaY", i] * myStates["sigmaY", i+1])
        J = 1/2
        Hint = J * Hint # change interaction strength

        H0 = 0
        for i in range(0,N+1):
                H0 = H0 + myStates["sigmaZ", i]
        Hges = Hint + H0
        return Hges

    elif interaction=="long-range":
        Hint=0
        for k in range(0,N+1):
            for i in range(0,k):
                    J = (1/(np.abs(i-k)))**1.22
                    Hint = Hint + J * (myStates["sigmaP", k] * myStates["sigmaM", i] + myStates["sigmaP", i] * myStates["sigmaM", k])
        H0 = 0
        for i in range(0,N+1):
            H0 = H0 + myStates["sigmaZ", i]
        Hges = Hint + H0
        return Hges



def timeEvo(dt, rho, Hint): #time evolution of an operator rho
	# The more effective way to make expontial of a matrix is using 'scipy.linalg.expm'
	U = -1.j * Hint * dt
	U = U.expm() 
	#U = la.expm(-1.j * Hint * dt)
	# In this case it's faster to define the unitary dagger than apply 'matrix.getH'(= dagger)
	Ud = 1.j * Hint * dt
	Ud = Ud.expm()
	#Ud = la.expm(1.j * Hint * dt)
	#return np.matmul(Ud, np.matmul(rho, U))
	rho = rho * Ud
	return U * rho


"""
def fidelity(rho, sigma):
	sqrt_rho = la.sqrtm(rho)
	argument = np.matmul(sqrt_rho, np.matmul(sigma, sqrt_rho))
	return (np.trace(la.sqrtm(argument)))**2
	# fidelity is usually defined for pure states, this definition uses density matrices
	# so that mixed states can be parsed as well.
	# https://doi.org/10.1016/0034-4877(76)90060-4
"""


def qubitfidelity(rho, sigma):# only for 2x2 matrices
	tr = np.trace(np.matmul(rho,sigma)) + 2 * np.sqrt(la.det(rho) * la.det(sigma))


def stateTransfer(dt, Hint, first, last):
	U = la.expm(-1.j * Hint * dt)
	Fdings = np.matmul(U,first)
	F = np.matmul(last, Fdings)
	return F
	

def deriv(dt, f):
	fdot = [(f[i+1]-f[i])/(dt) for i in range(len(f)-1)]
	return fdot


def integrate(N, psi0, tlist, solver, interaction):

    si = qeye(2)
    sx = sigmax()
    sy = sigmay()
    sz = sigmaz()

    sx_list = []
    sy_list = []
    sz_list = []

    for n in range(N+1):
        op_list = [si for m in range(N+1)]

        op_list[n] = sx
        sx_list.append(tensor(op_list))

        op_list[n] = sy
        sy_list.append(tensor(op_list))

        op_list[n] = sz
        sz_list.append(tensor(op_list))

    # construct the hamiltonian
    H = generate_states(N, interaction)

    # energy splitting terms

    # interaction terms
    

    # collapse operators
    c_op_list = []


    # evolve and calculate expectation values
    if solver == "me":
        result = mesolve(H, psi0, tlist, c_op_list, sz_list)
    elif solver == "mc":
        ntraj = 250 
        result = mcsolve(H, psi0, tlist, c_op_list, sz_list, ntraj)
    return result.expect