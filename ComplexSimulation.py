import numpy as np
import MonteCarlo
from datetime import datetime
from multiprocessing import Pool, cpu_count
import itertools
import multiprocessing as mp
import copy
import types

STATUS = {'Exit': 1, 'Died': 0}

def _run_simulation_worker(simulation):
    simulation.run_simulation()
    return simulation

def _parallel_simulations(simulations_list, num_processes=None):

    if num_processes is None:
        num_sim = len(simulations_list)
        num_processes = mp.cpu_count()
        if num_processes < num_sim:

            raise ValueError('Num of sim more than num of spu')
    
    simulations_copies = [copy.deepcopy(sim) for sim in simulations_list]

    with mp.Pool(processes=num_processes) as pool:

        results_sims = pool.map(_run_simulation_worker, simulations_copies)
    
    return results_sims

def _copy_func(f, name=None):
    return copy.deepcopy(f)

def save_spectrum(list_of_res, file_name):

    res_file_QE = open(f'{file_name}_QE.txt', 'w')
    res_file_p_intensity = open(f'{file_name}_p_intensity.txt', 'w')
    res_file_p_transport = open(f'{file_name}_p_transport.txt', 'w')
    res_file_emmitance = open(f'{file_name}_emmitance.txt', 'w')

    for dict_res in list_of_res:
        gamma = dict_res['gamma']
        QE = dict_res['QE']
        Emmitance = dict_res['Emmitance']
        N_el_exit = dict_res['N_el_exit']
        R = dict_res['R']

        res_file_QE.write(f'{gamma}\t{QE}\t{Emmitance}\t{N_el_exit}\t{R}\n')
        
        res_file_p_intensity.write(f'{gamma}\t')
        for indx, p_inten in enumerate(dict_res['p_intensity']):

            if indx != len(dict_res['p_intensity']) - 1:
                res_file_p_intensity.write(f'{p_inten}\t')
            else:
                res_file_p_intensity.write(f'{p_inten}\n')
        
        res_file_p_transport.write(f'{gamma}\t')
        for indx, p_trans in enumerate(dict_res['p_transports']):
            if indx != len(dict_res['p_transports']) - 1:
                res_file_p_transport.write(f'{p_trans}\t')
            else:
                res_file_p_transport.write(f'{p_trans}\n')

        res_file_emmitance.write(f'{gamma}\t')
        for indx, k_em in enumerate(dict_res['part_emmitance']):
            if indx != len(dict_res['part_emmitance']) - 1:
                res_file_emmitance.write(f'{k_em}\t')
            else:
                res_file_emmitance.write(f'{k_em}\n')


    res_file_QE.close()
    res_file_p_intensity.close()
    res_file_p_transport.close()
    res_file_emmitance.close()

def _save_exp_res(res: dict, exp_name: str):

    with open(f'res_{exp_name}.txt', 'w') as f:

        f.write(f'Result of simulation with name {exp_name}\n')
        current_datetime = datetime.now()
        f.write(f'Simulation perfomed at {str(current_datetime)}\n')

        for key in res:

            f.write(f'{key} = {res[key]}\n')

# only K2CsSb!!!!
def _calculate_emmitance(all_exit_electron, all_emmitance):
    if all_exit_electron != 0:
        return np.sqrt(0.117*all_emmitance/all_exit_electron)*np.sqrt(3.2/(9.1*2.99))*1e-2
    else:
        return 0

class Complex_Simulation:

    def __init__(self, gamma: float, N_sim: int, exp_name = '', parallel_flag = False):

        self.gamma = gamma
        self.simulation_mass = []
        self.N_sim = N_sim
        self.exp_name = exp_name
        self.parallel_flag = parallel_flag

        for i in range(N_sim):
            if parallel_flag:
                self.simulation_mass.append(MonteCarlo.Simulation(gamma, copy.deepcopy(exp_name) + f'_indx={int(i)}', i))
            else:
                self.simulation_mass.append(MonteCarlo.Simulation(gamma, copy.deepcopy(exp_name) + f'_indx={int(i)}'))
        
    def set_semiconductor(self, semiconductor_name, semiconductor):

        for sim in self.simulation_mass:
            semiconductor_temp = copy.deepcopy(semiconductor)
            semiconductor_name_temp = copy.deepcopy(semiconductor_name)
            sim.set_semiconductor(semiconductor_name_temp, semiconductor_temp)

    def set_calc_params(self, dt, N_subsim, N_el_in_subsim, N_iterations, kill_energy):

        self.N_subsim = N_subsim
        self.N_el_in_subsim = N_el_in_subsim
        self.N_el = len(self.simulation_mass)*self.N_subsim*self.N_el_in_subsim

        for sim in self.simulation_mass:

            sim.set_calc_params(dt, N_subsim, N_el_in_subsim, N_iterations, kill_energy)

    def set_sim_params(self, mass_dict_params: list, energy_DOS):

        self.mass_dict = mass_dict_params

        if len(self.mass_dict) != self.N_sim:

            raise ValueError('len(mass_dict_params) != len(self.simulation_mass)')

        for indx in range(self.N_sim):

            curr_sim = self.simulation_mass[indx]
            curr_sim.set_DOS(copy.deepcopy(energy_DOS), copy.deepcopy(mass_dict_params[indx]['coor_DOS']))
            curr_sim.set_geometry(copy.deepcopy(mass_dict_params[indx]['geometry']))

    def add_l_e_e_scattering(self, scat_func, delta_E):

        for sim in self.simulation_mass:

            sim.add_l_e_e_scattering(_copy_func(scat_func), delta_E)

    def run_simulation(self):

        if self.parallel_flag:
            self.simulation_mass = _parallel_simulations(self.simulation_mass)
        else:
            for sim in self.simulation_mass:

                sim.run_simulation()
    
    def get_results(self):

        QE = 0
        p_transports = []
        p_intensity = []
        emmitance_mass = []
        N_el_exit = 0
        k_quads = 0

        R = self.mass_dict[0]['R']

        for indx in range(self.N_sim):

            curr_sim = self.simulation_mass[indx]
            curr_dict_param = self.mass_dict[indx]
            p_exit = curr_sim.get_results()
            Abs = (1-curr_dict_param['R'])
            p = curr_dict_param['p']
            emmitance_mass.append(curr_sim.get_emittance())  

            QE += Abs*p*p_exit
            N_el_exit += p_exit*self.N_subsim*self.N_el_in_subsim

            p_transports.append(p_exit)
            p_intensity.append(curr_dict_param['p'])
            k_quads += curr_sim.get_emittance()

        emmitance = 1000*_calculate_emmitance(N_el_exit, k_quads)

        res = {'gamma': self.gamma,
              'QE': QE,
             'Emmitance': emmitance,
             'p_intensity': p_intensity,
             'p_transports': p_transports,
             'N_el_exit': N_el_exit,
             'R': R,
             'part_emmitance': emmitance_mass}

        _save_exp_res(res, self.exp_name)

        return res


