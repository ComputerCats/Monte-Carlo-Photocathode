import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_x_distr(name_pict, electron_gas):

    N_electrons = electron_gas.shape[0]

    for i in range(N_electrons):
        plt.scatter(electron_gas[i, 0], electron_gas[i, 2], color = 'blue')

    plt.grid()
    plt.xlabel('X')
    plt.ylabel('Y')

    plt.savefig(f'slide={name_pict}.png')

def plot_E_distr(name_pict, electron_gas):

    fig, ax = plt.subplots()

    ax.scatter(electron_gas[:, 2], electron_gas[:, -1], color = 'blue')

    ax.grid()
    ax.set_ylabel('E')
    ax.set_xlabel('Z')

    plt.savefig(f'E_slide={name_pict}.png')

def plot_energy_history(name_pict, history_energy, gamma, E_g):

    fig, ax = plt.subplots()

    ax.plot(history_energy[:, 0], history_energy[:, 1], color = 'red', label = f'gamma = {gamma}, E_g = {E_g}')
    ax.legend()
    ax.grid()
    ax.set_xlabel('Time, fs')
    ax.set_ylabel('Energy, eV')

    plt.savefig(f'{name_pict}.png')

def plot_coor_distr(name_pict, electron_gas):

    dz = 0.003

    min_z = np.min(electron_gas[:, 2])
    max_z = np.max(electron_gas[:, 2])

    N_iter = int((max_z - min_z)/dz) - 1

    fig, ax = plt.subplots()

    for i in range(N_iter):

        prev_down = min_z + i*dz
        curr_up = min_z + (i+1)*dz

        prev_mask = electron_gas[:, 2] > prev_down
        curr_mask = electron_gas[:, 2] < curr_up

        all_mask = np.logical_and(prev_mask, curr_mask)

        curr_value_plot = electron_gas[all_mask, 2].shape[0]

        ax.scatter(prev_down, curr_value_plot, color = 'blue')

    ax.grid()
    ax.set_ylabel('N electrons')
    ax.set_xlabel('Z')

    fig.savefig(f'{name_pict}.png')

def plot_initial_energy_distr(name_pict, electron_gas, all_energies):

    fig, ax = plt.subplots()
    
    N_electrons = all_energies.shape[0]

    result = np.zeros((N_electrons, 2))
    
    for i in range(N_electrons):

        energy = all_energies[i]

        N_electron_energy = electron_gas[electron_gas[:, -1] == energy, -1].shape[0]

        result[i, 0] = energy
        result[i, 1] = N_electron_energy

    ax.plot(result[:, 0], result[:, 1], color = 'blue')
    ax.grid()
    ax.set_ylabel('N electrons')
    ax.set_xlabel('E, ev')

    fig.savefig(f'{name_pict}.png')

def plot_coor_dos(file_name, coor_dos):

    fig, ax = plt.subplots()

    ax.plot(coor_dos[:, 2], coor_dos[:, 3])

    ax.grid()
    ax.set_xlabel('z, $\mu$')
    ax.set_ylabel('Coor DOS')

    fig.savefig(f'{file_name}')

def compare_with_exp(way_to, QE):

    exp_data = pd.read_csv(way_to, header = None, sep = '; ').to_numpy()

    fig, ax = plt.subplots()

    ax.plot(QE[:, 0], 100*QE[:, 1], label = 'Monte Carlo', color = 'red')
    ax.scatter(exp_data[:, 0], exp_data[:, 1], label = 'Experiment', color = 'blue')

    ax.grid()
    ax.set_xlabel('$\hbar\omega$')
    ax.set_ylabel('QE')
    ax.legend()

    fig.savefig('Comparing.png')

def compare_with_exp(way_to_exp, QE, way_to_val):

    exp_data = pd.read_csv(way_to_exp, header = None, sep = '; ').to_numpy()
    val_data = pd.read_csv(way_to_val, header = None, sep = '; ').to_numpy()

    fig, ax = plt.subplots()

    ax.plot(QE[:, 0], 100*QE[:, 1], label = 'Monte Carlo', color = 'red')
    ax.scatter(exp_data[:, 0], exp_data[:, 1], label = 'Experiment', color = 'blue')
    ax.scatter(val_data[:, 0], val_data[:, 1], label = 'Validation', color = 'black')

    ax.grid()
    ax.set_xlabel('$\hbar\omega$')
    ax.set_ylabel('QE')
    ax.legend()

    fig.savefig('Comparing.png')