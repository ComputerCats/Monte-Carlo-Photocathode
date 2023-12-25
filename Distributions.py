import numpy as np
import pandas as pd
from scipy import interpolate
import Geometry

def _make_electron_DOS(way_to, E_g, sep_file = '; '):

    #in this data assumed: E_f = 0

    electron_dos = pd.read_csv(way_to, sep = sep_file, header = None)
    electron_dos.iloc[:, 1] = electron_dos.iloc[:, 1] - np.min(electron_dos.iloc[:, 1]) #assume min DOS must be zero 
    energyes = electron_dos.iloc[:, 0].to_numpy() - E_g/2
    dos = electron_dos.iloc[:, 1].to_numpy()
    func = interpolate.interp1d(energyes, dos)

    return func

def _make_electron_DOS_Cs3Sb(way_to, E_g, sep_file = '; '): #Cs3Sb

    #in this data assumed: E_v = 0

    electron_dos = pd.read_csv(way_to, sep = sep_file, header = None)
    electron_dos.iloc[:, 1] = electron_dos.iloc[:, 1] - np.min(electron_dos.iloc[:, 1]) #assume min DOS must be zero 
    energyes = electron_dos.iloc[:, 0].to_numpy() - E_g
    dos = electron_dos.iloc[:, 1].to_numpy()
    func = interpolate.interp1d(energyes, dos)

    return func

def make_energy_DOS(way_to, E_g, gamma, delta_E, sep_file = '; '):

    DOS_func = _make_electron_DOS(way_to, E_g, sep_file)

    N_energyes = int((gamma - E_g)/delta_E)

    internal_integ_func = lambda E, gamma: DOS_func(E)*DOS_func(E-gamma)

    norm = Geometry.L_2_norm(internal_integ_func, gamma, E_g, gamma)

    eletrons_energyes = np.zeros((N_energyes, 2))

    for i in range(N_energyes):

        eletrons_energyes[i, 0] = i*delta_E
        eletrons_energyes[i, 1] = DOS_func(eletrons_energyes[i, 0])*DOS_func(eletrons_energyes[i, 0] - gamma)*delta_E

    eletrons_energyes[:, -1] = eletrons_energyes[:, -1]/norm

    delta_norm = (1 - np.cumsum(eletrons_energyes[:, -1], axis = 0)[-1])/N_energyes
    eletrons_energyes[:, -1] += delta_norm
    eletrons_energyes[:, 0] = eletrons_energyes[:, 0]

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
