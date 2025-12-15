import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

sys.path.append(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')
sys.path.append(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\mat_prop')

import MatProp
import Geometry
import MonteCarlo
import Distributions
from scipy import interpolate
import MyScatterings as scat
import Visualization

def test_make_DOS_func(way_to):

    res_func = Distributions._make_electron_DOS(way_to, 1.6, sep_file = '; ')

    energyes = np.linspace(0, 1.2, 100)
    DOS_1 = res_func(energyes)
    DOS_2 = res_func(energyes - 2.5)

    fig, ax = plt.subplots()

    ax.plot(energyes, DOS_1, label = 'DOS_0', color = 'blue')
    ax.plot(energyes, DOS_2, label = 'DOS_minus', color = 'red')
    ax.plot(energyes, DOS_1*DOS_2, label = 'DOS_reaction', color = 'black')
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

def test_make_energy_DOS_Cs3Sb():

    way_to_en_DOS = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\experiment\Cs3Sb\Cs3Sb_DOS.csv'

    E_g = 1.6 #band gap
    delta_E = 0.001

    energy_DOS = Distributions.make_energy_DOS(way_to_en_DOS, E_g, 2.3, delta_E)

    print(f'result_DOS = {energy_DOS}')
    print(f'norm must be equal 1, curr sum = {np.cumsum(energy_DOS[:, -1], axis = 0)[-1]}')

    plt.plot(energy_DOS[:, 0], energy_DOS[:, 1])
    plt.grid()
    plt.xlabel('Energy above cond band')
    plt.ylabel('DOS')
    plt.savefig('energy_dos.png')

def test_energy_and_coor_hist():

    way_to_en_DOS = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\mat_prop\data\K2CsSb\Ettema_k2cssb_dos.csv'

    N_electrons = 20000
    N_iterations = 20000
    E_g = 1.1 #band gap
    E_g_dos = 1.2
    E_a = 0.8 #electron afinity
    delta_E = 0.027 #ev, phonon energy
    delta_E_DOS = 0.01 #ev
    effective_mass = 0.117 #m/m_e
    kill_energy = E_a/1.5 #ev
    dt = 0.1
    semiconductor = scat.Semiconductor(E_a, E_g, effective_mass)

    z0 = 0.512

    K2CsSb = MatProp.K2CsSb()

    coor_DOSes = K2CsSb.get_coor_and_gamma_from_fdtd_for_plate()

    h = 0.040

    fig, ax = plt.subplots()

    N_iters = 100000

    gamma_cur = 2.32
    task = MonteCarlo.Simulation(gamma_cur, exp_name = f'test_energy')

    coor_DOS = K2CsSb.get_halfspace_coor_distr(gamma_cur, N_points = 400, z1 = z0 + h, z0 = z0)
    energy_DOS = Distributions.make_energy_DOS(way_to_en_DOS, E_g_dos, gamma_cur, delta_E_DOS)

    task.set_semiconductor('K2CsSb', semiconductor)
    task.set_calc_params(dt, N_electrons, N_iterations, kill_energy)
    task.set_DOS(energy_DOS, coor_DOS)
    task.add_l_e_e_scattering(scat.l_e_e, -E_g)
    task.add_l_e_e_scattering(scat.l_phonon, -delta_E)
    geom = Geometry.PlateGeom(np.array([0, 0, 0.512 + h]), np.array([0, 0, 1]), np.array([0, 0, 0.512]), np.array([0, 0, -1]))
    task.set_geometry(geom)

    result_dos = np.zeros((N_iters, 2))

    for i in range(N_iters):

        print(f'progres = {i/N_iters}')
        single_electron = task._initial_process_single_electron()
        result_dos[i, 0] = single_electron.get_E()
        result_dos[i, 1] = single_electron.get_coor()[-1]
    
    print('finish')

    nbins = 200

    fig, ax = plt.subplots(2, figsize = (10, 10))

    ax[0].hist(result_dos[:, 0], bins=nbins, edgecolor='black')
    ax[1].hist(result_dos[:, 1], bins=nbins, edgecolor='blue')

    ax[0].grid()
    ax[1].grid()

    ax[0].set_title('Energy distr')
    ax[1].set_title('Coor distr')

    fig.savefig('DOSes')

#test_make_DOS_func(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\experiment\K2CsSb\K2CsSb_DOS.csv')
#test_make_DOS_func(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\experiment\Cs3Sb\Cs3Sb_DOS.csv')
#test_make_energy_DOS_Cs3Sb()
test_energy_and_coor_hist()