import MyScatterings as scat
import numpy as np
import matplotlib.pyplot as plt

def test_scat():

    E_interval = np.linspace(0.35, 1.3, 100)

    fig, ax = plt.subplots()

    for energy in E_interval:

        curr_tau_1 = scat.tau_POP_plus(energy)
        curr_tau_2 = scat.tau_POP_minus(energy)

        ax.scatter(energy, curr_tau_1, color = 'blue')
        ax.scatter(energy, curr_tau_2, color = 'red')

    ax.grid()
    ax.set_title('Blue - plus, red - minus')

    fig.savefig('distr/Tau_and_energy.png')

test_scat()
