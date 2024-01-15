import numpy as np
import pandas as pd
from scipy import interpolate
import Geometry
import Distributions
import ElTransport
import Validation as val
import electron
import ElectronExit

EXIT_STATUS = Geometry.STATUS

def make_calc_log(i, exited_electrons, initial_electrons):

    print(f'Calculation progress: {i/initial_electrons}')
    print(f'Curr Yield: {exited_electrons/initial_electrons}')


class Simulation:

    #way_to - way to save
    #way_from - way to library

    def __init__(self, way_from, way_to, gamma):

        self.way_from = way_from
        self.way_to = way_to
        self.gamma = gamma

        self.scatterings_tau = []
        self.scatterings_E = []

        self._init_log_mass()

    def get_val_hitory_states(self):

        return np.array(self.history_electron_states)

    def _init_log_mass(self):

        self.mass_str_log = []

    def _add_str_to_log(self, other_str):

        self.mass_str_log.append(other_str)

    def set_semiconductor(self, semiconductor):

        self.semiconductor = semiconductor

    #dt fs
    def set_calc_params(self, dt, N, N_iterations, kill_energy):

        self.N_iterations = N_iterations
        self.initial_N_electrons = N
        self.dt = dt
        self.kill_energy = kill_energy

    #tau fs
    def add_scattering(self, tau, delta_E):

        self.scatterings_tau.append(tau)
        self.scatterings_E.append(delta_E)

    def set_geometry(self, geometry):

        self.geometry = geometry
    
    def initial_process(self, energy_DOS, coor_DOS):

        #electron columns = [x, y, z, phi (0, 2pi), psi (0, pi), E]
        #energy_DOS = columns: (energy, propability)
        #coor_DOS = columns: (x, y, z, propability)

        numbers_of_position = range(0, coor_DOS.shape[0])

        #prepare positions
        indx_pos = np.random.choice(numbers_of_position, p=coor_DOS[:, -1].reshape(1, -1)[0])

        #set energyes, dir and coor
        coor = coor_DOS[indx_pos,:3]
        energy = np.random.choice(energy_DOS[:,0], p=energy_DOS[:,1])
        direction = ElTransport.make_initial_dir()

        #Visualization.plot_coor_distr(r'distr\initial_electon_coor_distr', self.electron_gas)
        #Visualization.plot_initial_energy_distr(r'distr\initial_electon_energy_distr', self.electron_gas, self.energyes_DOS[:,0])

        single_electron = electron.Electrons(coor[0], coor[1], coor[2], direction[0], direction[1], energy)

        return single_electron

    def run_simulation(self, energy_DOS, coor_DOS):

        self.exit_electron = 0

        for i in range(self.initial_N_electrons):

            make_calc_log(i, self.exit_electron, self.initial_N_electrons)

            #Visualization.plot_x_distr(str(i), self.electron_gas)
            #Visualization.plot_E_distr(str(i), self.electron_gas)
            self._run_new_iteration(energy_DOS, coor_DOS)

    def kill_low_energy_electrons(self, single_electron, kill_energy):

        return single_electron.get_E > kill_energy

    def _run_new_iteration(self, energy_DOS, coor_DOS):

        single_electron = self.initial_process(energy_DOS, coor_DOS)

        for i in range(self.N_iterations):

            if self.kill_low_energy_electrons(single_electron, self.kill_energy):
                
                break

            ElTransport.transport_process(single_electron, self.tau_mass, self.E_mass, self.dt)
            
            if ElectronExit.exit_process(self.geom, single_electron):

                self.exit_electron += 1
                break


    def get_results(self):

        return self.exit_electron/self.initial_N_electrons

    def save_params_of_simulation(self):

        pass
        #log_file = open('', 'w')