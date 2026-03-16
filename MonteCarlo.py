import numpy as np
import ElTransport
import electron
import ElectronExit
import log
import copy

def make_console_log(i, exited_electrons, N_el, N_subsim):

    print(f'Calculation progress: {100*round((i+1)/N_subsim, 3)} %')
    print(f'Curr Yield: {round(exited_electrons/(N_el)*100, 1)} %')

class Simulation:

    #way_to - way to save
    #way_from - way to library

    def __init__(self, gamma, exp_name = '', sim_indx = None):

        self.exp_name = exp_name

        self.gamma = gamma

        self._init_scat_mass()

        self.log_exp = log.MyLog(self.exp_name)

        self.sim_indx = sim_indx

        self.is_save_xy_coor = False

    def set_save_xy_coor(self):

        self.is_save_xy_coor = True
        self.xy_mass = []

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

    def _init_N_el_in_subsims(self, N_el, N_el_in_subsim):

        self.mass_N_el = [N_el_in_subsim for i in range(N_el // N_el_in_subsim)]
        self.mass_N_el.append(N_el % N_el_in_subsim)

    #dt fs
    def set_calc_params(self, dt, N_el, N_iterations, kill_energy, N_el_in_subsim):

        self.N_iterations = N_iterations
        self.kill_energy = kill_energy
        self.dt = dt
        self.n_electrons_in_subsim = N_el_in_subsim
        self.barrier_flag = True
        self.N_el = N_el

        self._init_N_el_in_subsims(N_el, N_el_in_subsim)

    def add_l_e_e_scattering(self, l_e_e, delta_E):

        self.scatterings_l_e_e.append(l_e_e)
        self.scatterings_E_l_e_e.append(delta_E)

    def get_sim_indx(self):

        return self.sim_indx

    def set_geometry(self, geometry):

        self.geometry = geometry

    def _initial_process_electron(self, N_el):

        #electron columns = [x, y, z, vx, vy, vz, {0 or 1}]

        electrons = electron.Electrons(N_el)
        electrons.set_electron_properties(self.semiconductor.get_effective_mass())

        numbers_of_position = range(0, self.coor_DOS.shape[0])
        for indx_el in range(N_el):
            indx_pos = np.random.choice(numbers_of_position, p=self.coor_DOS[:, -1].reshape(1, -1)[0])

            coor = self.coor_DOS[indx_pos,:3]
            energy = np.random.choice(self.energy_DOS[:,0], size = (N_el, 1), p=self.energy_DOS[:,1])

            electrons.set_xcoor(np.array([coor[0], coor[1], coor[2]]), indx_el)
        
        veloicity = ElTransport.make_initial_veloicity(electrons, energy)

        return electrons

    def run_simulation(self):

        self.exit_electron = 0
        self.emmitance = 0
        self.N_subsim = len(self.mass_N_el)

        for i in range(self.N_subsim):
           
            self._run_new_subsim(self.mass_N_el[i])
            if self.sim_indx is None:
                make_console_log(i, self.exit_electron, self.N_el, self.N_subsim)

        self._end_experiment()

    def _run_new_subsim(self, N_el):
        electrons = self._initial_process_electron(N_el)

        for i in range(self.N_iterations): #time

            electrons.kill_low_energy_electron(self.kill_energy)
            if electrons.is_end(): break

            ElTransport.transport_process(electrons, self.dt, self.scatterings_l_e_e, self.scatterings_E_l_e_e)
            res_iter = ElectronExit.exit_process(self.geometry, electrons, self.semiconductor, self.barrier_flag, self.is_save_xy_coor)

            self.exit_electron += res_iter['N_exit']
            self.emmitance += res_iter['Emmitance']

            if self.is_save_xy_coor:

                new_mass = np.zeros((3, 1))
                new_mass = {'coor': res_iter['xy_mass'], 'time': self.dt*i}

                self.xy_mass.append(new_mass)

            if self.sim_indx is None:
                if i == self.N_iterations -1:
                    print('All time')

    def save_xy_mass(self, file_name = 'xy_distribution'):

        xy_distribution_file = open(f'{file_name}.txt', 'w')

        xy_distribution_file.write('x\ty\ttime\n')

        if self.is_save_xy_coor:

            N_iter = len(self.xy_mass)

            for i in range(N_iter):

                curr_iter_mass = self.xy_mass[i]
                coor = self.xy_mass[i]['coor']
                time = self.xy_mass[i]['time']
                N_coor_in_iter = len(coor)

                for j in range(N_coor_in_iter):
                    
                    xy_distribution_file.write(f'{coor[j][0]}\t{coor[j][1]}\t{time}\n')
        else:

            raise ValueError('XY mass is empty')

        xy_distribution_file.close()
                
    def _add_params_to_log(self):

        self.log_exp._add_str_to_log('dt', f'{self.dt}')
        self.log_exp._add_str_to_log('N_iterations', f'{self.N_iterations}')
        self.log_exp._add_str_to_log('n_electrons_in_subsim', f'{self.mass_N_el}')
        self.log_exp._add_str_to_log('N_subsim', f'{len(self.mass_N_el)}')
        self.log_exp._add_str_to_log('kill_energy', f'{self.kill_energy}')

        self.log_exp._add_str_to_log('geometry', f'{self.geometry.get_name()}, params = {self.geometry.get_params()}')

        self.log_exp._add_str_to_log('semiconductor', f'{self.semiconductor_name}, E_a = {self.semiconductor.get_E_a()}, E_g = {self.semiconductor.get_E_g()}')

        for indx, delta_E in enumerate(self.scatterings_E_l_e_e):

            self.log_exp._add_str_to_log('scattering ', f'delta E = {delta_E} l (0.5 ev) = {self.scatterings_l_e_e[indx](0.5)}\n')

        self.log_exp._add_str_to_log('QE_=_', f'{self.get_results()}')
        self.log_exp._add_str_to_log('Emmitance_=_', f'{self.emmitance}')

    def _end_experiment(self):
        
        self._add_params_to_log()
        self.log_exp.save_log()

    def get_emittance(self):

        return self.emmitance
    
    def get_results(self):

        return {'N photoelectrons': self.exit_electron, 'N_el': self.N_el, 'P transport': self.exit_electron/self.N_el}

    
        