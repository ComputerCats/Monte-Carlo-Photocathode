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
import Validation as val
import Visualization

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

def plot_spectrum(way_to):

    way_to = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\spectrums\Cs3Sb'
    way_from = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\spectrums\Cs3Sb'
    way_to_en_DOS = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\experiment\Cs3Sb\Cs3Sb_DOS.csv'

    N_electrons = 10000
    N_iterations = 10000
    E_g = 1.6 #band gap
    E_a = 0.3 #electron afinity
    delta_E = 0.027 #ev, phono energy
    delta_E_DOS = 0.002 #ev
    effective_mass = 0.12 #m/m_e
    dt = 20 #fs
    kill_energy = E_a/2 #ev

    fig, ax = plt.subplots()

    N_energies = 1

    energyes = np.linspace(1.95, 2.6, N_energies)

    R_func = get_R_func()

    results = np.zeros((N_energies, 2))

    for i in range(N_energies):

        gamma_cur = energyes[i]
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

        validation.plot_energy_DOS(r'\Validation\energy', energy_DOS)
        validation.plot_coor_DOS(r'\Validation\coor', coor_DOS)
        validation.plot_Z_electron_distr(r'Validation\distr', coor_DOS)
        validation.plot_initial_energy_distr(r'Validation\distr', energy_DOS[:, 0])

        task.run_simulation()

        refl = R_func(gamma_cur)

        results[i, 0] = gamma_cur
        results[i, 1] = (1-refl)*task.get_results()

    ax.plot(results[:, 0], results[:, 1], label = 'QE')

    ax.grid()
    ax.set_xlabel('$\hbar\omega$')
    ax.set_ylabel('$QE$')

    fig.savefig('result_spectr1.png')
    with open(f'result_spectr1.npy', 'wb') as f:
        np.save(f, results)

def plot_ready_results(file_name):

    with open(f'{file_name}.npy', 'rb') as f:
        result = np.load(f)

    Visualization.compare_with_exp(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\experiment\Cs3Sb\ExpCs3Sb.csv', result)

E_photon_interval = np.linspace(1.95, 3, 10)
#plot_spectrum(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')
plot_ready_results('result_spectr1')
