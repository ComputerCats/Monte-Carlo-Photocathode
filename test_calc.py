import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import Geometry

sys.path.append(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')

import electron


def test_init_phase_prostr(way_to):

    way_to = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode'

    N_electrons = 10000
    N_iterations = 10000
    E_g = 1.2
    E_a = 0.7
    gamma = 2.5
    delta_E = 0.01
    delta_E_DOS = 0.002
    l_E = 0.003

    task = electron.Electrons(way_to, gamma)
    task.set_calc_params(N_electrons, N_iterations)
    task.set_semi_params(E_g, E_a, delta_E_DOS)
    task.set_transport_params(delta_E, l_E)
    task.initial_phase_prostr()
    initial_electron_gas = np.copy(task.electron_gas)
    geom = Geometry.HalfspaceGeom(np.array([0, 0, 0.512]), np.array([0, 0, 1]))
    print('initial_gas = ' + str(initial_electron_gas))
    task.set_geometry(geom)
    task.run_simulation()
    print(f'result_distr = {task.electron_gas}')
    print(f'N_out = {task.exit_electron}')
    print(f'result = {task.get_results()}')

test_init_phase_prostr(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')