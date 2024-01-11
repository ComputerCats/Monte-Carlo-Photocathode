import numpy as np
import pandas as pd
from scipy import interpolate
import Geometry

def _make_electron_DOS(way_to, E_g, sep_file = '; '):

    #in this data assumed: E_f = 0

    electron_dos = pd.read_csv(way_to, sep = sep_file, header = None)
    electron_dos.iloc[:, 1] = electron_dos.iloc[:, 1] 
    energyes = electron_dos.iloc[:, 0].to_numpy() - E_g
    dos = np.abs(electron_dos.iloc[:, 1].to_numpy())
    func = interpolate.interp1d(energyes, dos)

    return func

def make_energy_DOS(way_to, E_g, gamma, delta_E, sep_file = '; '):

    DOS_func = _make_electron_DOS(way_to, E_g, sep_file)

    N_energyes = int((gamma - E_g)/delta_E)

    norm = 0

    eletrons_energyes = np.zeros((N_energyes, 2))

    for i in range(N_energyes):

        eletrons_energyes[i, 0] = i*delta_E
        eletrons_energyes[i, 1] = DOS_func(eletrons_energyes[i, 0])*DOS_func(eletrons_energyes[i, 0] - gamma)*delta_E
        norm += DOS_func(eletrons_energyes[i, 0])*DOS_func(eletrons_energyes[i, 0] - gamma)*delta_E

    eletrons_energyes[:, -1] = eletrons_energyes[:, -1]/norm

    delta_norm = (1 - np.cumsum(eletrons_energyes[:, -1], axis = 0)[-1])
    print(f'delta_norm = {delta_norm}')
    eletrons_energyes[-1, 1] += delta_norm

    return eletrons_energyes

def make_coordinate_DOS(way_to):

    coordinates = pd.read_csv(way_to, sep = '\t')
    net = coordinates.iloc[:, 1:4].to_numpy()
    field = coordinates.iloc[:, 4].to_numpy()

    intensity = field*field

    norm = np.cumsum(intensity, axis = 0)[-1]
    result_DOS = intensity.reshape(-1, 1)/norm
    result = np.concatenate((net, result_DOS), axis=1)

    return result
