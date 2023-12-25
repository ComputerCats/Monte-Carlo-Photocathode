import numpy as np

# life of electron: (x, y, z, psi, theta, E) -> phonon scattering (change angle, change coor, change energy) -> new iter


def make_new_coor(electron_gas, dt, effective_mass): 

    l_E = 0.001*np.sqrt(electron_gas[:, -1]*2*1.6/(effective_mass*9.1))*dt

    cos_theta_mass = np.cos(electron_gas[:, 4])
    sin_theta_mass = np.sin(electron_gas[:, 4])
    sin_psi_mass = np.sin(electron_gas[:, 3])
    cos_psi_mass = np.cos(electron_gas[:, 3])

    electron_gas[:, 0] += l_E*sin_theta_mass*cos_psi_mass
    electron_gas[:, 1] += l_E*sin_theta_mass*sin_psi_mass
    electron_gas[:, 2] += l_E*cos_theta_mass

    return electron_gas

def _make_p_mass(E, dt, tau):

    if tau(E) <= 0:
        
        raise 'Tau must be greater then 0'

    if tau(E) >= 1:
        
        raise 'dt/tau must be less then 1'

    p = dt/tau(E)

    return [p, 1-p]

# refactor
def make_new_energy(electron_gas, tau_mass, E_mass, dt):

    N_electrons = electron_gas[:, 0].shape[0]

    for i in range(N_electrons):

        for j in range(len(tau_mass)):

            p_mass = _make_p_mass(electron_gas[j, -1], dt, tau_mass[j])

            delta_E = np.random.choice([E_mass[j], 0], p = p_mass)

            electron_gas[i, -1] =  electron_gas[i, -1] + delta_E

    return electron_gas

def make_new_dir(electron_gas):

    N_electrons = electron_gas.shape[0]

    electron_gas[:, 3:5] = np.random.rand(N_electrons, 2)

    electron_gas[:, 3] = 2*np.pi*electron_gas[:, 3]
    electron_gas[:, 4] = np.pi*electron_gas[:, 4]

    return electron_gas