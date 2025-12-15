import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import Geometry

sys.path.append(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')

import electron

def test_init_phase_prostr(way_to):

    way_to = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode'

    N_electrons = 10
    N_iterations = 10
    E_g = 1.2
    E_a = 0.7
    gamma = 2.5
    delta_E = 0.027
    delta_E_DOS = 0.005

    task = electron.Electrons(way_to, N_electrons, N_iterations, E_g, E_a, gamma, delta_E, delta_E_DOS)
    task.initial_phase_prostr()
    initial_electron_gas = np.copy(task.electron_gas)
    geom = Geometry.HalfspaceGeom(np.array([0, 0, 0.512]), np.array([0, 0, 1]))
    print('initial_gas = ' + str(initial_electron_gas))
    task.set_geometry(geom)
    task._run_new_iteration()
    new_electron_gas = task.electron_gas
    print(f'new_iter = {new_electron_gas}')
    print(f'N_out = {task.exit_electron}')

test_init_phase_prostr(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')