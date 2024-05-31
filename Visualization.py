import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate


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

def compare_with_exp(way_to, QE, title):

    exp_data = pd.read_csv(way_to, header = None, sep = '; ').to_numpy()

    fig, ax = plt.subplots()

    ax.plot(QE[:, 0], 100*QE[:, 1], label = 'Monte Carlo', color = 'red')
    ax.scatter(exp_data[:, 0], exp_data[:, 1], label = 'Experiment', color = 'blue')

    ax.grid()
    ax.set_xlabel('$\hbar\omega$, eV')
    ax.set_ylabel('QE, %')
    ax.set_xlim(left = 2.0, right = 2.4)
    ax.set_ylim(top = 14)
    ax.set_title(title)
    ax.legend()

    fig.savefig('Comparing.png')

def compare_with_exp_and_val(way_to_exp, QE, way_to_val):

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

def compare_with_another_result(QE, another_res, this_res_name ='', another_res_name='', title=''):

    fig, ax = plt.subplots()

    ax.plot(QE[:, 0], 100*QE[:, 1], label = this_res_name, color = 'blue')
    ax.scatter(QE[:, 0], 100*QE[:, 1], color = 'blue')

    ax.plot(another_res[:, 0], 100*another_res[:, 1], label = another_res_name, color = 'red')
    ax.scatter(another_res[:, 0], 100*another_res[:, 1], color = 'red')

    ax.grid()
    ax.set_xlabel('$\hbar\omega$, eV')
    ax.set_ylabel('QE, %')
    ax.legend()
    ax.set_title(title)
    ax.set_xlim(left = 1.8, right = 2.4)
    ax.set_ylim(top = 16)

    fig.savefig('Comparing.png')

def compare_with_exp_another_result(way_to_exp, QE, another_res, this_res_name ='', another_res_name='', title=''):

    fig, ax = plt.subplots()

    exp_data = pd.read_csv(way_to_exp, header = None, sep = '; ').to_numpy()

    ax.plot(QE[:, 0], 100*QE[:, 1], label = this_res_name, color = 'red')
    ax.scatter(QE[:, 0], 100*QE[:, 1], color = 'red')

    ax.plot(another_res[:, 0], 100*another_res[:, 1], label = another_res_name, color = 'blue')
    ax.scatter(another_res[:, 0], 100*another_res[:, 1], color = 'blue')

    ax.plot(exp_data[:, 0], exp_data[:, 1], label = 'Experiment', color = 'k')
    ax.scatter(exp_data[:, 0], exp_data[:, 1], color = 'k')

    ax.grid()
    ax.set_xlabel('$\hbar\omega$, eV')
    ax.set_ylabel('QE, %')
    ax.legend()
    ax.set_title(title)
    ax.set_xlim(left = 2.0, right = 2.4)
    ax.set_ylim(top = 14)

    fig.savefig('Comparing.png')

def compare_three_res(res1, res2, res3, res1_name ='',res2__name='', res3_name ='', title=''):

    fig, ax = plt.subplots()

    ax.plot(res1[:, 0], 100*res1[:, 1], label = res1_name, color = 'red')
    ax.scatter(res1[:, 0], 100*res1[:, 1], color = 'red')

    ax.plot(res2[:, 0], 100*res2[:, 1], label = res2__name, color = 'blue')
    ax.scatter(res2[:, 0], 100*res2[:, 1], color = 'blue')

    ax.plot(res3[:, 0], 100*res3[:, 1], label = res3_name, color = 'k')
    ax.scatter(res3[:, 0], 100*res3[:, 1], color = 'k')

    ax.grid()
    ax.set_xlabel('$\hbar\omega$, eV')
    ax.set_ylabel('QE, %')
    ax.legend()
    ax.set_title(title)
    ax.set_xlim(left = 2.0, right = 2.4)
    ax.set_ylim(top = 16)

    fig.savefig('Comparing.png')

def plot_p(res, res_name = '', title = ''):

    fig, ax = plt.subplots()

    ax.plot(1000*res[:, 0], 100*res[:, 1], label = res_name, color = 'red')
    ax.scatter(1000*res[:, 0], 100*res[:, 1], color = 'red')

    ax.grid()
    ax.set_xlabel('h, nm')
    ax.set_ylabel('Propability, %')
    ax.set_xlim(left = 14)
    ax.legend()
    ax.set_title(title)

    fig.savefig('Comparing.png')

def plot_error(dict_QE, etalon):

    fig, ax = plt.subplots()

    for dt in dict_QE:

        error = np.max(abs(etalon[:, 1] - dict_QE[dt][:, 1])/etalon[:, 1])
        ax.scatter(float(dt), error, color = 'red')

    ax.set_xlabel('dt, fs')
    ax.set_ylabel('error, %')
    ax.grid()

    fig.savefig('error.png')

def plot_ratio(QE1, QE2, xlabel = '', ylabel = '', title = ''):

    func_QE1 = interpolate.interp1d(QE1[:,0], QE1[:,1])
    func_QE2 = interpolate.interp1d(QE2[:,0], QE2[:,1])

    energies = np.linspace(2, 2.4, 100)

    fig, ax = plt.subplots()

    ratio = func_QE1(energies)/func_QE2(energies)

    ax.plot(energies, ratio, color = 'red')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid()
    ax.set_title(title)

    fig.savefig('Ratio.png')
