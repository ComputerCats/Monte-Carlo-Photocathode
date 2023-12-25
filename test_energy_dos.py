import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import Distributions

sys.path.append(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')

import electron

def test_make_DOS_func(way_to):

    res_func = Distributions._make_electron_DOS(way_to, 1.6, sep_file = '; ')

    energyes = np.linspace(-2, 2, 100)
    DOS = res_func(energyes)

    fig, ax = plt.subplots()

    ax.plot(energyes, DOS, label = 'DOS', color = 'blue')
    ax.grid()
    ax.set_xlabel('E, eV')
    ax.set_ylabel('DOS')
    ax.legend()

    fig.savefig('distr\Energy_DOS.png')

def test_make_energy_DOS(way_to):

    #test energy distribution
    energy_DOS = electron._make_energy_DOS(way_to, 1.6, 2, 0.002)

    print(f'result_DOS = {energy_DOS}')
    print(f'norm must be equal 1, curr sum = {np.cumsum(energy_DOS[:, -1], axis = 0)[-1]}')

    plt.plot(energy_DOS[:, 0], energy_DOS[:, 1])
    plt.grid()
    plt.xlabel('Energy above cond band')
    plt.ylabel('DOS')
    plt.savefig('energy_dos.png')

#test_make_energy_DOS(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\experiment\Cs3Sb\Cs3Sb_DOS.csv')
test_make_DOS_func(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\experiment\Cs3Sb\Cs3Sb_DOS.csv')