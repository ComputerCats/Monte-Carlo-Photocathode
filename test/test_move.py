import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')

import Geometry
import MonteCarlo
import electron as el
import Distributions
from scipy import interpolate
import MyScatterings as scat
import Visualization
import Validation as val

def get_coor_Cs3Sb(gamma):
    way_to = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\experiment\Cs3Sb\alpha_Cs3Sb.csv'
    data = pd.read_csv(way_to, sep = '; ').to_numpy()

    func = interpolate.interp1d(data[:, 0], data[:, 1])

    N_points = 400

    curr_alpha = func(gamma)

    coor_dos = np.zeros((N_points, 4))

    dx = (0.912-0.512)/N_points

    for i in range(N_points):

        coor_dos[i, 3] = np.exp(-curr_alpha*dx*i)
        coor_dos[i, 2] = i*dx + 0.512

    norm = np.cumsum(coor_dos[:, 3], axis = 0)[-1]
    coor_dos[:, 3] = coor_dos[:, 3]/norm

    Visualization.plot_coor_dos('Coor_distribution.png', coor_dos)

    return coor_dos

def get_curr_gamma(way_to_coor_DOS, file_name):
    way_to = f'{way_to_coor_DOS}\\{file_name}'
    coordinates = pd.read_csv(way_to, sep = '\t')

    wave_length = coordinates.iloc[0, 0]
    gamma = 1/wave_length*2*np.pi*3*1/(16)

    return gamma

def get_R_func():
    
    func = lambda x: 0.2

    return func

def test_plot_electron_history():

    way_to = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\spectrums\Cs3Sb'
    way_from = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\spectrums\Cs3Sb'
    way_to_en_DOS = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\experiment\Cs3Sb\Cs3Sb_DOS.csv'

    N_electrons = 1
    N_iterations = 10000
    E_g = 1.6 #band gap
    E_a = 0.3 #electron afinity
    delta_E = 0.027 #ev, phono energy
    delta_E_DOS = 0.001 #ev
    effective_mass = 0.12 #m/m_e
    dt = 4 #fs
    kill_energy = E_a/2 #ev

    N_iter = 10

    energy = 2.5

    for i in range(N_iter):

        gamma_cur = energy
        print(f'gamma_cur = {gamma_cur}')
        coor_DOS = get_coor_Cs3Sb(gamma_cur)
        energy_DOS = Distributions.make_energy_DOS(way_to_en_DOS, E_g, gamma_cur, delta_E)
        
        electron = el.Electrons()
        electron.set_semi_params(E_g, E_a, delta_E_DOS)
        electron.set_electron_propities(effective_mass)

        task = MonteCarlo.Simulation(way_from, way_to, gamma_cur)
        task.set_electron(electron)
        task.set_calc_params(dt, N_electrons, N_iterations, kill_energy)
        task.initial_phase_prostr(energy_DOS, coor_DOS)

        task.add_scattering(scat.tau_POP_plus, delta_E)
        task.add_scattering(scat.tau_POP_minus, -delta_E)

        geom = Geometry.HalfspaceGeom(np.array([0, 0, 0.512]), np.array([0, 0, 1]))
        task.set_geometry(geom)

        validation = val.ValidateSim(task)
        validation.plot_energy_DOS(energy_DOS)
        validation.plot_initial_energy_distr(energy_DOS[:, 0])

        task.set_validation()

        task.run_simulation()

        validation.plot_electron_states(f'electron_history_indx={i}', task.get_val_hitory_states())

test_plot_electron_history()
