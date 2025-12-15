import MyScatterings as scat
import numpy as np
import matplotlib.pyplot as plt

def test_scat():

    E_interval = np.linspace(0.2, 1.5, 100)

    fig, ax = plt.subplots()

    for energy in E_interval:

        curr_tau_1 = scat.tau_POP_plus(energy)
        curr_tau_2 = scat.tau_POP_minus(energy)

        ax.scatter(energy, curr_tau_1, color = 'blue')
        ax.scatter(energy, curr_tau_2, color = 'red')

    ax.grid()
    ax.set_title('Blue - plus, red - minus')

    fig.savefig('distr/Tau_and_energy.png')

def test_scat_l_e():

    fig, ax = plt.subplots()

    N = 100

    E_interval = np.linspace(1.3, 2, N)

    data = np.zeros((N, 2))

    for indx in range(N):

        data[indx, 0] = E_interval[indx]
        data[indx, 1] = scat.l_e_e(E_interval[indx])

    ax.scatter(E_interval, data[:, 1], color = 'blue')

    ax.grid()
    ax.set_title('l_e_e')
    ax.set_xlabel('$\mu$')
    ax.set_ylabel('$l_e$')

    fig.savefig('distr/l_e_e.png')

test_scat()
test_scat_l_e()
