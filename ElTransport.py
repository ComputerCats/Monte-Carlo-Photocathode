import numpy as np

# life of electron: (x, y, z, psi, theta, E) -> phonon scattering (change angle, change coor, change energy) -> new iter


def _make_new_coor(single_electron, dt): 

    l_E = single_electron.get_veloicity()*dt

    curr_dir = single_electron.get_dir()

    cos_theta_mass = np.cos(curr_dir[1])
    sin_theta_mass = np.sin(curr_dir[1])
    sin_psi_mass = np.sin(curr_dir[0])
    cos_psi_mass = np.cos(curr_dir[0])

    dx = l_E*sin_theta_mass*cos_psi_mass
    dy = l_E*sin_theta_mass*sin_psi_mass
    dz = l_E*cos_theta_mass

    single_electron.add_coor(np.array([dx, dy, dz]))

def _make_p_mass(E, dt, tau):

    p = 1 - np.exp(-dt/tau(E))

    if tau(E) <= 0:
        
        raise ValueError('Tau must be greater then 0')

    if p >= 1:
        
        raise ValueError('dt/tau must be less then 1')

    return [p, 1-p]

def _make_scatterings(single_electron, tau_mass, E_mass, dt):

    N_tau = len(tau_mass)

    for j in range(N_tau):

        p_mass = _make_p_mass(single_electron.get_E(), dt, tau_mass[j])

        is_scattering = np.random.choice([True, False], p = p_mass)

        if is_scattering:

            single_electron.add_energy(E_mass[j])

            _make_new_dir(single_electron)

def _make_new_dir(single_electron):

    new_psi = 2*np.pi*np.random.rand()
    new_theta = np.pi*np.random.rand()

    single_electron.set_dir(np.array([new_psi, new_theta]))

def make_initial_dir():

    new_psi = 2*np.pi*np.random.rand()
    new_theta = np.pi*np.random.rand()

    return np.array([new_psi, new_theta])

def transport_process(single_electron, tau_mass, E_mass, dt):

    _make_new_coor(single_electron, dt)
    _make_scatterings(single_electron, tau_mass, E_mass, dt)

