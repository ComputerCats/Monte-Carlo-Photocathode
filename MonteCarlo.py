import numpy as np
import ElTransport
import electron
import ElectronExit
import log

EXIT_STATUS= ElectronExit.EXIT_PROCESS_STATUS

def make_console_log(i, exited_electrons, initial_electrons):

    print(f'Calculation progress: {100*round(i/initial_electrons, 3)} %')
    print(f'Curr Yield: {round(exited_electrons/initial_electrons*100, 1)} %')

class Simulation:

    #way_to - way to save
    #way_from - way to library

    def __init__(self, gamma, exp_name = ''):

        self.exp_name = exp_name

        self.gamma = gamma

        self._init_scat_mass()

        self.log_exp = log.MyLog(self.exp_name)

    def _init_scat_mass(self):

        self.scatterings_l_e_e = []
        self.scatterings_E_l_e_e = []

    def set_semiconductor(self, semiconductor_name, semiconductor):

        self.semiconductor = semiconductor
        self.semiconductor_name = semiconductor_name
        
    def set_DOS(self, energy_DOS, coor_DOS):

        #energy_DOS = columns: (energy, propability)
        #coor_DOS = columns: (x, y, z, propability)

        self.energy_DOS = energy_DOS
        self.coor_DOS = coor_DOS

    #dt fs
    def set_calc_params(self, dt, N, N_iterations, kill_energy):

        self.N_iterations = N_iterations
        self.initial_N_electrons = N
        self.kill_energy = kill_energy
        self.dt = dt

    def add_l_e_e_scattering(self, l_e_e, delta_E):

        self.scatterings_l_e_e.append(l_e_e)
        self.scatterings_E_l_e_e.append(delta_E)

    def set_geometry(self, geometry):

        self.geometry = geometry

    def _initial_process_single_electron(self):

        #electron columns = [x, y, z, phi (0, 2pi), psi (0, pi), E]

        numbers_of_position = range(0, self.coor_DOS.shape[0])

        indx_pos = np.random.choice(numbers_of_position, p=self.coor_DOS[:, -1].reshape(1, -1)[0])

        coor = self.coor_DOS[indx_pos,:3]
        energy = np.random.choice(self.energy_DOS[:,0], p=self.energy_DOS[:,1])

        single_electron = electron.Electrons(coor[0], coor[1], coor[2], 0, 0, 0)
        single_electron.set_electron_properties(self.semiconductor.get_effective_mass())

        veloicity = ElTransport.make_initial_dir(single_electron, energy)
        single_electron.set_veloicity(veloicity)

        return single_electron

    def run_simulation(self):

        self.exit_electron = 0
        self.emmitance = 0

        for i in range(self.initial_N_electrons):

            make_console_log(i, self.exit_electron, self.initial_N_electrons)

            self._run_new_iteration()

        self._end_experiment()

    def _get_this_electron_emitance(self, single_electron):

        cos_out = self.geometry.get_cos_angle(single_electron)

        return single_electron.get_E()*(1-cos_out*cos_out)

    def _calculate_mean(self, ):

        if self.exit_electron != 0:
            return np.sqrt(self.semiconductor.get_effective_mass()*self.emmitance/self.exit_electron)*np.sqrt(3.2/(9.1*3))*1e-2
        else:
            return 0

    def _run_new_iteration(self):
        single_electron = self._initial_process_single_electron()

        for i in range(self.N_iterations):

            ElTransport.transport_process(single_electron, self.dt, self.scatterings_l_e_e, self.scatterings_E_l_e_e)
            electron_status = ElectronExit.exit_process(self.geometry, single_electron, self.semiconductor, self.kill_energy)

            if electron_status == EXIT_STATUS['Out']:
                self.emmitance += self._get_this_electron_emitance(single_electron)
                self.exit_electron += 1
                break

            if electron_status == EXIT_STATUS['Died']:

                break

    def _add_params_to_log(self):

        self.log_exp._add_str_to_log('dt', f'{self.dt}')
        self.log_exp._add_str_to_log('N_iterations', f'{self.N_iterations}')
        self.log_exp._add_str_to_log('initial_N_electrons', f'{self.initial_N_electrons}')
        self.log_exp._add_str_to_log('kill_energy', f'{self.kill_energy}')

        self.log_exp._add_str_to_log('geometry', f'{self.geometry.get_name()}, params = {self.geometry.get_params()}')

        self.log_exp._add_str_to_log('semiconductor', f'{self.semiconductor_name}, E_a = {self.semiconductor.get_E_a()}, E_g = {self.semiconductor.get_E_g()}')

        for indx, delta_E in enumerate(self.scatterings_E_l_e_e):

            self.log_exp._add_str_to_log('scattering ', f'delta E = {delta_E} l (0.5 ev) = {self.scatterings_l_e_e[indx](0.5)}\n')

        self.log_exp._add_str_to_log('QE_=_', f'{self.get_results()}')
        self.log_exp._add_str_to_log('Emmitance_=_', f'{self.final_emmitance}')

    def _end_experiment(self):
        
        self.final_emmitance = self._calculate_mean()
        self._add_params_to_log()
        self.log_exp.save_log()

    def get_emittance(self):

        return self.final_emmitance
    
    def get_results(self):

        return self.exit_electron/self.initial_N_electrons

    
        