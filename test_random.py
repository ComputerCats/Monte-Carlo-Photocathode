import ElTransport as transp
import numpy as np
import matplotlib.pyplot as plt
import MyScatterings as myscat

def _check_res(mat_real, mat_exp):

    return abs(1-mat_exp/mat_real)


def write_res(log_file, E, dt, is_good, sign_of_tau):

    log_file.write(f'{is_good}\t{E}\t{dt}\t{sign_of_tau}\n')

def make_log_file(header_info):

    log_file = open(r'test\test_random\log.txt', 'w')

    log_file.write(header_info)

    return log_file

def test_choice():

    T_all = 1e5
    dt_mass = [2, 3, 4, 5, 6, 7, 8, 9, 10]

    tau_plus = myscat.tau_POP_plus
    tau_minus = myscat.tau_POP_minus

    E_mass = np.linspace(0.5, 1.9, 10)

    log_file = make_log_file(f'T_all={T_all}, dt_mass={dt_mass}\nis_good\tE\tdt\tsign_of_tau\n')

    for dt in dt_mass:

        print(f'dt={dt}\n')

        for E in E_mass:

            N_exp = int(T_all/dt)

            mat_plus = 0
            mat_minus = 0

            p_mass_plus = transp._make_p_mass(E, dt, tau_plus)
            p_mass_minus = transp._make_p_mass(E, dt, tau_minus)

            for i in range(N_exp):

                mat_plus += np.random.choice([1, 0], p = p_mass_plus)
                mat_minus += np.random.choice([1, 0], p = p_mass_minus)

            mat_real_plus = N_exp*(1 - np.exp(-dt/tau_plus(E)))
            mat_real_minus = N_exp*(1 - np.exp(-dt/tau_minus(E)))

            write_res(log_file, round(E, 2), dt, round(_check_res(mat_real_plus, mat_plus), 2), '+')
            write_res(log_file, round(E, 2), dt, round(_check_res(mat_real_minus, mat_minus), 2), '-')

    log_file.close()
    

test_choice()        


                

