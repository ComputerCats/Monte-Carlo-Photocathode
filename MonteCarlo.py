import numpy as np
import pandas as pd
from scipy import interpolate
import Geometry
import Distributions
import ElTransport

EXIT_STATUS = Geometry.STATUS

def p_exit(E, E_a, cos_angle):

    E_exit = E*cos_angle*cos_angle

    if E <= 0:

        return 0

    if (E_exit <= E_a) or (np.sqrt(E_a/E) < cos_angle):

        return 0

    result = 4*np.sqrt(E_exit*(E_exit-E_a))/(np.sqrt(E_exit-E_a)+np.sqrt(E_exit))**2
    return result


def make_calc_log(i, N_iterations, Curr_n_electrons, exited_electrons, initial_electrons):

    print(f'Calculation progress: {i/N_iterations}')
    print(f'Count of electron in volume: {Curr_n_electrons}')
    print(f'Curr Yield: {exited_electrons/initial_electrons}')


class Simulation:

    #way_to - way to save
    #way_from - way to library

    def __init__(self, way_from, way_to, gamma):

        self.way_from = way_from
        self.way_to = way_to
        self.gamma = gamma

        self.scatterings_p = []
        self.scatterings_E = []

        self._init_log_mass()

    def _init_log_mass(self):

        self.mass_str_log = []

    def _add_str_to_log(self, other_str):

        self.mass_str_log.append(other_str)

    def set_electron(self, electron):

        self.electron = electron

        self._add_str_to_log(self.electron.get_log())

    def _init_distributions(self, way_to_en_DOS, way_to_coor_DOS):

        self.energyes_DOS = Distributions._make_energy_DOS(way_to_en_DOS, self.electron.E_g, self.gamma, self.electron.delta_E_DOS)
        self.pos_DOS = Distributions._make_coordinate_DOS(way_to_coor_DOS)

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

    def _prepare_scat(self):

        all_scatt_prop = 0

        for x in self.scatterings_p:

            all_scatt_prop += x

        if x > 1:

            raise Exception('Propability of scattering > 1. Check your dt and tau\'s')

        prop_zero_scat = 1-all_scatt_prop

        #zero scat
        self.scatterings_p.append(prop_zero_scat)
        self.scatterings_E.append(0)


    def set_geometry(self, geometry):

        self.geometry = geometry
    
    def initial_phase_prostr(self, energy_DOS, coor_DOS, data_type = '1d'):

        #electron gas columns = [x, y, z, phi (0, 2pi), psi (0, pi), E]
        #energy_DOS = columns: (energy, propability)
        #coor_DOS = columns: (x, y, z, propability)

        self.electron_gas = np.zeros((self.initial_N_electrons, 6))
        self.electron_gas = ElTransport.make_new_dir(self.electron_gas)

        numbers_of_position = range(0, coor_DOS.shape[0])

        for i in range(self.initial_N_electrons):

            #set positions
            indx_pos = np.random.choice(numbers_of_position, p=coor_DOS[:, -1].reshape(1, -1)[0])

            #set energyes
            self.electron_gas[i, :3] = coor_DOS[indx_pos,:3]
            self.electron_gas[i, 5] = np.random.choice(energy_DOS[:,0], p=energy_DOS[:,1])

        #Visualization.plot_coor_distr(r'distr\initial_electon_coor_distr', self.electron_gas)
        #Visualization.plot_initial_energy_distr(r'distr\initial_electon_energy_distr', self.electron_gas, self.energyes_DOS[:,0])

    def run_simulation(self):

        self._prepare_scat()

        self.exit_electron = 0

        for i in range(self.N_iterations):

            Curr_n_electrons = self.electron_gas.shape[0]

            make_calc_log(i, self.N_iterations, Curr_n_electrons, self.exit_electron, self.initial_N_electrons)

            if Curr_n_electrons == 0:
                break

            #Visualization.plot_x_distr(str(i), self.electron_gas)
            #Visualization.plot_E_distr(str(i), self.electron_gas)
            self._run_new_iteration()

    def exit_and_kill(self, electron_gas):

        N_electrons = self.electron_gas.shape[0]

        kill_list = []

        for i in range(N_electrons):

            status = self.geometry.get_status(self.electron_gas[i, :3])

            if status == EXIT_STATUS['Exit']:
                prop_exit = p_exit(self.electron_gas[i, -1], self.electron.E_a, self.geometry.get_cos_angle([self.electron_gas[i, 3], self.electron_gas[i, 4]]))
                is_exit = np.random.choice([False, True], p = [1-prop_exit, prop_exit])

                if is_exit:

                    kill_list.append(i)

                    #important string!!!
                    self.exit_electron += 1
                else:
                    self.electron_gas[i, 3:5] = self.geometry.reflect(self.electron_gas[i, 3:5])

        self.electron_gas = np.delete(self.electron_gas, kill_list, axis = 0)

    def kill_low_energy_electrons(self, kill_energy):

        self.electron_gas = self.electron_gas[self.electron_gas[:, -1] > kill_energy, :]

    def _run_new_iteration(self):

        self.kill_low_energy_electrons(self.kill_energy)
        self.electron_gas = ElTransport.make_new_coor(self.electron_gas, self.dt, self.electron.effective_mass)
        self.electron_gas = ElTransport.make_new_dir(self.electron_gas)
        self.electron_gas = ElTransport.make_new_energy(self.electron_gas, self.scatterings_tau, self.scatterings_E)
        self.electron_gas = self.exit_and_kill(self.electron_gas)

    def get_results(self):

        return self.exit_electron/self.initial_N_electrons

    def save_params_of_simulation(self):

        pass
        #log_file = open('', 'w')