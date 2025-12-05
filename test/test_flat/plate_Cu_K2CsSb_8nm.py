import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

sys.path.append(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')
sys.path.append(r'G:\kintech\Diplom\mat_prop')

import MatProp
import Geometry
import MonteCarlo
import Distributions
from scipy import interpolate
import MyScatterings as scat
import Visualization

STATUS = {'Exit': 1, 'Died': 0}

def l_phonon_plus(E):
    return 0.012

def l_phonon_minus(E):
    return 0.004

def plot_spectrum():

    way_to_en_DOS = r'G:\kintech\Diplom\mat_prop\data\K2CsSb\K2CsSb_DOS.csv'

    N_subsims = 100
    N_el_in_subsim = 100
    N_iterations = 20000
    E_g = 1.1 #band gap
    E_a = 0.7 #electron afinity
    delta_E = 0.027 #ev, phonon energy
    delta_E_DOS = 0.01 #ev
    effective_mass = 0.117 #m/m_e
    kill_energy = 5*delta_E #ev
    dt = 0.1
    semiconductor = scat.Semiconductor(E_a, E_g, effective_mass)

    z0 = 0.512
    h = 0.008

    K2CsSb = MatProp.K2CsSb_plate()
    K2CsSb.set_fdtd_way(r'G:\kintech\Diplom\mat_prop\data\K2CsSb\plate_fdtd\8nmCu', 'Cuprum_epsilon_fdtd.txt')

    cases = K2CsSb.get_coor_and_gamma_from_fdtd_for_plate(0.512)

    results = np.zeros((len(cases), 6))

    geom = Geometry.ConvexShape()

    plane_cathode = Geometry.Plane(np.array([0, 0, 0.512-h]), np.array([0, 0, -1]), 'True')
    plane_down = Geometry.Plane(np.array([0, 0, 0.512]), np.array([0, 0, 1]), 'plane_down')

    geom.add_plane(plane_cathode, STATUS['Exit'])
    geom.add_plane(plane_down, STATUS['Died'])

    for i, gamma_cur in enumerate(cases):
        #if i != 0: continue
        print(f'gamma_cur = {gamma_cur}')
        
        coor_DOS = cases[gamma_cur]['Distr_box']
        energy_DOS = Distributions.make_energy_DOS(way_to_en_DOS, E_g, gamma_cur, delta_E_DOS)

        task = MonteCarlo.Simulation(gamma_cur, f'gamma_cur_=_{gamma_cur}')
        task.set_semiconductor('K2CsSb', semiconductor)
        task.set_calc_params(dt, N_subsims, N_el_in_subsim, N_iterations, kill_energy)
        task.set_DOS(energy_DOS, coor_DOS)

        task.add_l_e_e_scattering(l_phonon_plus, delta_E)
        task.add_l_e_e_scattering(l_phonon_minus, -delta_E)

        task.set_geometry(geom)

        task.run_simulation()

        results[i, 0] = gamma_cur
        results[i, 1] = (1-cases[gamma_cur]['R'])*cases[gamma_cur]['p']*task.get_results()
        results[i, 2] = cases[gamma_cur]['p']
        results[i, 3] = (1-cases[gamma_cur]['R'])
        results[i, 4] = task.get_results()
        results[i, 5] = task.get_emittance()

        print('finish')


    with open(f'result_spectr1.npy', 'wb') as f:
        np.save(f, results)

def plot_ready_results(file_name):

    with open(f'{file_name}.npy', 'rb') as f:
        result = np.load(f)
    with open(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\spectrums\K2CsSb\plate_new\24nm\Cu\result_spectr1.npy', 'rb') as f:
        another_result = np.load(f)

    Visualization.compare_with_another_result(result, another_result, this_res_name ='8 nm', another_res_name='24 nm')

def plot_this_res(file_name):

    with open(f'{file_name}.npy', 'rb') as f:
        result = np.load(f)

    fig, ax = plt.subplots()

    ax.plot(result[:, 0], 100*result[:, 1], color = 'red')
    ax.scatter(result[:, 0], 100*result[:, 1], color = 'red')

    ax.set_xlabel('$\hbar\omega$, eV')
    ax.set_ylabel('QE, %')
    ax.grid()

    fig.savefig('this_res')

plot_spectrum()
plot_ready_results('result_spectr1')
plot_this_res('result_spectr1')
