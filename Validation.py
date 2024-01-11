import Distributions
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Visualization

class ValidateSim:
    
    def __init__(self, simulation):

        self.simulation = simulation

    @staticmethod
    def plot_energy_DOS(E_DOS):

        fig, ax = plt.subplots()

        ax.plot(E_DOS[:,0], 100*E_DOS[:,1], color = 'blue')
        ax.grid()
        ax.set_xlabel('E, eV')
        ax.set_ylabel('DOS')

        fig.savefig('energy_DOS.png')

    def plot_electron_DOS(self, way_to):

        DOS_func = Distributions._make_electron_DOS(way_to, self.simulation.electron.E_g)
        E = np.linspace(-1.6, 2, 100)

        fig, ax = plt.subplots()

        ax.plot(E, 100*DOS_func(E), color = 'blue')
        ax.grid()
        ax.set_xlabel('E above cond band, eV')
        ax.set_ylabel('DOS, %')

        fig.savefig('electron_DOS.png')

    @staticmethod
    def plot_coor_DOS(coor_DOS):

        fig, ax = plt.subplots()

        ax.plot(coor_DOS[:,2], 100*coor_DOS[:,3], color = 'blue')
        ax.grid()
        ax.set_xlabel('$\Delta z$, $\mu$')
        ax.set_ylabel('DOS')

        fig.savefig('coor_DOS.png')

    @staticmethod
    def plot_electron_states(pict_name, states):

        fig, ax = plt.subplots(figsize = (10, 10))

        N_states = states.shape[0]

        for indx_state in range(N_states):

            ax.scatter(states[indx_state, 0], states[indx_state, 2], color = 'red')
            if indx_state > 0:
                if states[indx_state-1, -1] != states[indx_state, -1]:
                    ax.annotate(f'{indx_state}, E = {round(states[indx_state, -1], 3)}', (states[indx_state, 0], states[indx_state, 2]))
            else:
                ax.annotate(f'{indx_state}, E = {round(states[indx_state, -1], 3)}', (states[indx_state, 0], states[indx_state, 2]))

        ax.plot(states[:, 0], states[:, 2], color = 'blue')
        ax.grid()
        ax.set_xlabel('$X$, $\mu$')
        ax.set_ylabel('$Z$, $\mu$')

        fig.savefig(f'{pict_name}.png')

    def plot_Z_electron_distr(self):

        Visualization.plot_coor_distr('Z_electron_distr.png', self.simulation.electron_gas)

    def plot_initial_energy_distr(self, all_energies):

        Visualization.plot_initial_energy_distr('energy_electron_distr.png', self.simulation.electron_gas, all_energies)
