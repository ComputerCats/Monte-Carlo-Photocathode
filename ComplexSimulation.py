import numpy as np
import ElTransport
import electron
import ElectronExit
import log
import MonteCarlo
from datetime import datetime

def _save_exp_res(res: dict, exp_name: str):

    with open(f'res_{exp_name}.txt', 'w') as f:

        f.write(f'Result of simulation with name {exp_name}\n')
        current_datetime = datetime.now()
        f.write(f'Simulation perfomed at {str(current_datetime)}\n')

        for key in res:

            f.write(f'{key} = {self.mass_log[key]}\n')

# only K2CsSb!!!!
def _calculate_emmitance_gamma(all_emmitance, all_exit_electron):
    if all_exit_electron != 0:
        return np.sqrt(0.117*all_emmitance/all_exit_electron)*np.sqrt(3.2/(9.1*2.99))*1e-2
    else:
        return 0
    
def _calculate_emmitance(N_exit, k_quads):
    
    all_exit_electron = 0
    summ_k_quads = 0
        
    for k in k_quads:
        
        summ_k_quads += k
        
    return calculate_emmitance_gamma(summ_k_quads, all_exit_electron)

class Complex_Simulation:

    def __init__(self, gamma: float, N_sim: int, exp_name = ''):

        self.gamma = gamma
        self.simulation_mass = []
        self.N_sim = N_sim

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

        if len(mass_dict_params) != self.N_sim:

            raise ValueError('len(mass_dict_params) != len(self.simulation_mass)')

        for indx in range(len(mass_dict_params)):

            curr_sim = self.simulation_mass[indx]

            curr_sim.set_semiconductor(mass_dict_params['semiconductor_name'], mass_dict_params['semiconductor'])
            curr_sim.set_DOS(energy_DOS, mass_dict_params['coor_DOS'])

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

        for indx in range(self.N_sim):

            curr_sim = self.simulation_mass[indx]
            curr_dict_param = self.mass_dict[indx]['Coor_distr']
            p_exit = curr_sim.get_results()
            Abs = (1-curr_dict_param['Coor_distr']['R'])
            p = curr_dict_param['Coor_distr']['p']
            emmitance_mass.append(self.simulation_mass[indx].get_emittance())  

            QE += Abs*p*p_exit
            N_el_exit += p_exit*self.N_el
            p_transports.append(p_exit)
            p_intensity.append(curr_dict_param['Coor_distr']['p']))
            k_quads += self.simulation_mass[indx].get_emittance()
            emmitans_mass.append(self.simulation_mass[indx].get_emittance())

        emmitance = _calculate_emmitance(N_el_exit, k_quads)

        res = ['gamma': self.gamma,
              'QE': QE,
             'Emmitance': emmitance,
             'p_intensity': p_intensity,
             'p_transports': p_transports,
             'N_el_exit': N_el_exit,
             'R': R,
             'part_emmitance': emmitance_mass
             ]


