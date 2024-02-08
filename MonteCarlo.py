import numpy as np
import pandas as pd
import Geometry
import ElTransport
import Validation as val
import electron
import ElectronExit
import Visualization

EXIT_STATUS = Geometry.STATUS

def make_calc_log(i, exited_electrons, initial_electrons):

    print(f'Calculation progress: {round(i/initial_electrons, 3)}')
    print(f'Curr Yield: {round(exited_electrons/initial_electrons*100, 1)} %')

class Simulation:

    #way_to - way to save
    #way_from - way to library

    def __init__(self, gamma):

        self.gamma = gamma

        self._init_scat_mass()
        self._init_log_mass()

    def get_val_hitory_states(self):

        return np.array(self.history_electron_states)

    def _init_log_mass(self):

        self.mass_str_log = []

    def _init_scat_mass(self):

        self.scatterings_tau = []
        self.scatterings_E = []
        self.scatterings_l_e_e = []
        self.scatterings_E_l_e_e = []

    def _add_str_to_log(self, other_str):

        self.mass_str_log.append(other_str)

    def set_semiconductor(self, semiconductor):

        self.semiconductor = semiconductor

    def set_DOS(self, energy_DOS, coor_DOS):

        self.energy_DOS = energy_DOS
        self.coor_DOS = coor_DOS

    #dt fs
    def set_calc_params(self, l_E, E_loss, N, N_iterations, kill_energy):

        self.N_iterations = N_iterations
        self.initial_N_electrons = N
        self.l_E = l_E
        self.kill_energy = kill_energy
        self.E_loss = E_loss

    #tau fs
    def add_scattering(self, tau, delta_E):

        self.scatterings_tau.append(tau)
        self.scatterings_E.append(delta_E)

    def add_l_e_e_scattering(self, l_e_e, delta_E):

        self.scatterings_l_e_e.append(l_e_e)
        self.scatterings_E_l_e_e.append(delta_E)

    def set_geometry(self, geometry):

        self.geometry = geometry
    
    def initial_process(self):

        #electron columns = [x, y, z, phi (0, 2pi), psi (0, pi), E]
        #energy_DOS = columns: (energy, propability)
        #coor_DOS = columns: (x, y, z, propability)

        numbers_of_position = range(0, self.coor_DOS.shape[0])

        #prepare positions
        indx_pos = np.random.choice(numbers_of_position, p=self.coor_DOS[:, -1].reshape(1, -1)[0])

        #set energyes, dir and coor
        coor = self.coor_DOS[indx_pos,:3]
        energy = np.random.choice(self.energy_DOS[:,0], p=self.energy_DOS[:,1])
        direction = ElTransport.make_initial_dir()

        single_electron = electron.Electrons(coor[0], coor[1], coor[2], direction[0], direction[1], energy)
        single_electron.set_electron_propities(self.semiconductor.get_effective_mass())

        return single_electron

    def run_simulation(self):

        self.exit_electron = 0

        for i in range(self.initial_N_electrons):

            make_calc_log(i, self.exit_electron, self.initial_N_electrons)

            self._run_new_iteration()

    def kill_low_energy_electrons(self, single_electron):

        if single_electron.get_E() < self.kill_energy:

            return True

        else:

            return False

    def _run_new_iteration(self):

        single_electron = self.initial_process()

        if self.kill_low_energy_electrons(single_electron): return

        for i in range(self.N_iterations):

            ElTransport.transport_process(single_electron, self.E_loss, self.l_E, self.scatterings_l_e_e, self.scatterings_E_l_e_e)

            if self.kill_low_energy_electrons(single_electron):
                
                break

            is_out = ElectronExit.exit_process(self.geometry, single_electron, self.semiconductor)
            
            if is_out:

                self.exit_electron += 1
                break

    def get_results(self):

        return self.exit_electron/self.initial_N_electrons

    def save_params_of_simulation(self):

        pass
        #log_file = open('', 'w')