import numpy as np
import MonteCarlo
from datetime import datetime


def save_spectrum(spectrum, name):

    with open(f'spectrum_{name}.txt', 'w') as f:

        f.write(f'gamma\tQE\tEmmitance\tR\n')

        for res in spectrum:

            gamma = res['gamma']
            QE = res['QE']
            Emmitance = res['Emmitance']
            R = res['R']

            f.write(f'{gamma}\t{QE}\t{Emmitance}\t{R}\n')

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

    def __init__(self, gamma: float, N_sim: int, exp_name = ''):

        self.gamma = gamma
        self.simulation_mass = []
        self.N_sim = N_sim
        self.exp_name = exp_name

        for i in range(N_sim):

            self.simulation_mass.append(MonteCarlo.Simulation(gamma, exp_name + f'_indx={int(i)}'))
        
    def set_semiconductor(self, semiconductor_name, semiconductor):

        for sim in self.simulation_mass:

            sim.set_semiconductor(semiconductor_name, semiconductor)

    def set_calc_params(self, dt, N, N_iterations, kill_energy):

        self.N_el = N

        for sim in self.simulation_mass:

            sim.set_calc_params(dt, N, N_iterations, kill_energy)

    def set_sim_params(self, mass_dict_params: list, energy_DOS):

        self.mass_dict = mass_dict_params

        if len(self.mass_dict) != self.N_sim:

            raise ValueError('len(mass_dict_params) != len(self.simulation_mass)')

        for indx in range(self.N_sim):

            curr_sim = self.simulation_mass[indx]
            curr_sim.set_DOS(energy_DOS, mass_dict_params[indx]['coor_DOS'])
            curr_sim.set_geometry(mass_dict_params[indx]['geometry'])

    def add_l_e_e_scattering(self, scat_func, delta_E):

        for sim in self.simulation_mass:

            sim.add_l_e_e_scattering(scat_func, delta_E)

    def run_simulation(self):

        for sim in self.simulation_mass:

            sim.run_simulation()

    def get_emittance(self):

        return self.final_emmitance
    
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
            N_el_exit += p_exit*self.N_el

            p_transports.append(p_exit)
            p_intensity.append(curr_dict_param['p'])
            k_quads += curr_sim.get_emittance()

        emmitance = _calculate_emmitance(N_el_exit, k_quads)

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


