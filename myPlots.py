from generateStates import *
from timeEvo import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import partialtrace
from qutip import *

# for using tex formatting and font in plots
#"""
plt.rcParams.update({"text.usetex": True, })
mpl.rcParams['text.latex.preamble'] = [
    r'\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{lmodern}\inputencoding{utf8}\usepackage{amsmath}\usepackage{amssymb}\usepackage{dsfont}\usepackage{mathtools}\usepackage{physics}']
mpl.rcParams['font.family'] = ['serif']


def plotZ(rho, save, title, t):  # i can expand this for x and y in the future
    Eges = np.zeros(len(t))
    energies = {}
    for i in range(N-1):
        energyZ = []
        for dt in t:
            rhoT = timeEvo(dt, rho, Hint)
            EoftZ = energy(rhoT.ptrace(i,i+1), sigmaz())
            energyZ.append(EoftZ)
        energies["energyZ", i] = energyZ
        for j in range(0, 200):
            Eges[j] += np.real(energyZ[j])

    # I chose hbar = 1 to avoid overflow errors. thats why it is time with units s/(Js) = 1/J.
    for i in range(N+1):
        zlabel = "$Tr[\sigma^z_{" + str(i) + "}\\rho]$"
        plt.plot(t, np.real(energies["energyZ", i]), label=zlabel)
        plt.title(title)
        plt.xlabel("time [1/J]")
        plt.ylabel("$Tr[\sigma^z\\rho]$")
        plt.legend()
        if save == True:
            plt.show(block=False)
            plt.savefig(title+".png")
        else:
            plt.show(block=False)
    plt.pause(10)


def plotFidelity(rho, t, title, saveas):
    #pt = partialtrace.PartialTrace(N)
    fidelityList = []
    rho_0 = rho.ptrace(0).unit()
    for dt in t:
        rho_X = timeEvo(dt, rho, Hint)
        rho_N = rho_X.ptrace(N).unit()
        fidelityList.append(fidelity(rho_0, rho_N))
    fidelityList = np.real(fidelityList)
    #print(fidelityList)
    fig, ax = plt.subplots()
    ax.plot(t, fidelityList, label="F(t)")#, '+')
    ax.set_title(title)
    ax.set_xlabel("time")
    ax.set_ylabel("Fidelity F(t)")
    ax.axvline(np.pi,color='grey')
    ax.axhline(1,color='grey')
    #plt.ylim(-0.2,1.2)
    plt.savefig(str(saveas) + ".png")
    plt.show(block=True)
    #plt.pause(5)


