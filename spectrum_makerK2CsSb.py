import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import Geometry
import MonteCarlo
import electron as el
import Distributions
from scipy import interpolate

sys.path.append(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')

import Visualization

def make_file_coor_name(indx):

    file_name_template = 'V_detector_f00'

    if indx < 10:

        file_name = f'{file_name_template}0{round(indx, 0)}.d'

    else:

        file_name = f'{file_name_template}{round(indx, 0)}.d'

    return file_name

def get_coor_DOS(way_to_coor_DOS, file_name):
    way_to = f'{way_to_coor_DOS}\\{file_name}'
    result = Distributions.make_coordinate_DOS(way_to)
    
    return result

def get_coor_Cs3Sb(gamma, N_points):
    way_to = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\experiment\Cs3Sb\alpha_Cs3Sb.csv'
    data = pd.read_scv(way_to, sep = '; ').to_numpy()

    func = interpolate.interp1d(data[:, 0], data[:, 1])

    curr_alpha = func(gamma)

    coor_dos = np.zeros((N_points, 4))

    coor_dos[:, 2] = (np.linspace(0.512, 0.912, 400)).reshape(1, -1)

    for i in range(400):

        coor_dos[i, 4] = np.exp(-curr_alpha*(coor_dos[i, 2] - 0.512))

    norm = np.cumsum(coor_dos[:, 4], axis = 0)[-1]
    coor_dos[:, 4] = coor_dos[:, 4]/norm

    return coor_dos

def get_curr_gamma(way_to_coor_DOS, file_name):
    way_to = f'{way_to_coor_DOS}\\{file_name}'
    coordinates = pd.read_csv(way_to, sep = '\t')

    wave_length = coordinates.iloc[0, 0]
    gamma = 1/wave_length*2*np.pi*3*1/(16)

    return gamma

def get_R_func(way_to_coor_DOS):

    flux_data = pd.read_csv(f'{way_to_coor_DOS}\\flux.d', sep = '\t').to_numpy()

    flux_data[:, 0] = 1/flux_data[:, 0]*2*np.pi*3*1/(16)
    reflection = flux_data[:, 1]/flux_data[:, 2]

    func = interpolate.interp1d(flux_data[:, 0], reflection)

    #func = lambda x: 0.2

    return func

def plot_spectrum(way_to):

    way_to = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode'
    way_from = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode'
    way_to_en_DOS = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\electronDOS.csv'
    way_to_coor_DOS = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\fdtd\halfspace'

    N_electrons = 80000
    N_iterations = 10000
    E_g = 1.1 #band gap
    E_a = 0.65 #electron afinity
    delta_E = 0.022 #ev, phono energy
    delta_E_DOS = 0.002 #ev
    tau = 26.6 #fs
    effective_mass = 0.12 #m/m_e
    dt = 5 #fs
    kill_energy = E_a #ev

    fig, ax = plt.subplots()

    N_energies = 16

    R_func = get_R_func(way_to_coor_DOS)
    
    results = np.zeros((N_energies, 2))

    for i in range(N_energies):

        file_coor_name = make_file_coor_name(i)
        gamma_cur = get_curr_gamma(way_to_coor_DOS, file_coor_name)
        print(f'gamma_cur = {gamma_cur}')
        coor_DOS = get_coor_DOS(way_to_coor_DOS, file_coor_name)
        energy_DOS = Distributions.make_energy_DOS(way_to_en_DOS, E_g, gamma_cur, delta_E)
        
        electron = el.Electrons()
        electron.set_semi_params(E_g, E_a, delta_E_DOS)
        electron.set_transport_params(delta_E, tau)
        electron.set_electron_propities(effective_mass)

        task = MonteCarlo.Simulation(way_from, way_to, gamma_cur)
        task.set_electron(electron)
        task.set_calc_params(dt, N_electrons, N_iterations, kill_energy)
        task.initial_phase_prostr(energy_DOS, coor_DOS)
        task.add_scatering()

        geom = Geometry.HalfspaceGeom(np.array([0, 0, 0.512]), np.array([0, 0, 1]))
        task.set_geometry(geom)
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

    Visualization.compare_with_exp(result)

E_photon_interval = np.linspace(1.95, 3, 10)
plot_spectrum(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')
plot_ready_results('result_spectr1')