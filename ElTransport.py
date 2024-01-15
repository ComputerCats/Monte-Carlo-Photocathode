import numpy as np

# life of electron: (x, y, z, psi, theta, E) -> phonon scattering (change angle, change coor, change energy) -> new iter

def get_electron_veloicity(energy, effective_mass):

    return 0.001*np.sqrt(energy*2*1.6/(effective_mass*9.1))

def make_new_coor(electron_gas, dt, effective_mass): 

    l_E = get_electron_veloicity(electron_gas[:, -1], effective_mass)*dt

    cos_theta_mass = np.cos(electron_gas[:, 4])
    sin_theta_mass = np.sin(electron_gas[:, 4])
    sin_psi_mass = np.sin(electron_gas[:, 3])
    cos_psi_mass = np.cos(electron_gas[:, 3])

    electron_gas[:, 0] += l_E*sin_theta_mass*cos_psi_mass
    electron_gas[:, 1] += l_E*sin_theta_mass*sin_psi_mass
    electron_gas[:, 2] += l_E*cos_theta_mass

    return electron_gas

def _make_p_mass(E, dt, tau):

    p = 1 - np.exp(-dt/tau(E))

    if tau(E) <= 0:
        
        raise ValueError('Tau must be greater then 0')

    if p >= 1:
        
        raise ValueError('dt/tau must be less then 1')

    return [p, 1-p]

def make_scatterings(electron_gas, tau_mass, E_mass, dt):

    N_electrons = electron_gas[:, 0].shape[0]
    N_tau = len(tau_mass)

    for i in range(N_electrons):

        for j in range(N_tau):

            p_mass = _make_p_mass(electron_gas[i, -1], dt, tau_mass[j])

            is_scattering = np.random.choice([True, False], p = p_mass)

            if is_scattering:

                electron_gas[i, -1] =  electron_gas[i, -1] + E_mass[j]

                electron_gas[i, :] = _make_new_dir(electron_gas[i, :])

    return electron_gas

def _make_new_dir(electron_gas):

    electron_gas[3] = 2*np.pi*np.random.rand()
    electron_gas[4] = np.pi*np.random.rand()

    return electron_gas

def initial_dir(electron_gas):

    N_electrons = electron_gas.shape[0]

    electron_gas[:, 3:5] = np.random.rand(N_electrons, 2)

    electron_gas[:, 3] = 2*np.pi*electron_gas[:, 3]
    electron_gas[:, 4] = np.pi*electron_gas[:, 4]

    return electron_gas