import numpy as np

# life of electron: (x, y, z, psi, theta, E) -> phonon scattering (change angle, change coor, change energy) -> new iter

def _make_new_coor(single_electron, dt): 

    l_E = 0.003#single_electron.get_veloicity()*dt

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

def _make_p_mass_l_e_e(single_electron, l_e_e, dt):

    p = 1 - np.exp(-single_electron.get_veloicity()*dt/l_e_e(single_electron.get_E()))

    if p >= 1:
        
        raise ValueError('dt/tau must be less then 1')

    return [p, 1-p]

def _make_scatterings(single_electron, tau_mass, E_mass):
    
    do_new_dir = False
    N_tau = len(tau_mass)

    for j in range(N_tau):

        single_electron.add_energy(E_mass[j])

        do_new_dir = True

    if do_new_dir:

        _make_new_dir(single_electron)

def _make_new_dir(single_electron):

    new_psi = 2*np.pi*np.random.rand()
    new_theta = np.pi*np.random.rand()

    single_electron.set_dir(np.array([new_psi, new_theta]))

def make_initial_dir():

    new_psi = 2*np.pi*np.random.rand()
    new_theta = np.pi*np.random.rand()

    return np.array([new_psi, new_theta])

def reflcation_process(geom, single_electron):

    new_coor = geom.get_new_point_after_reflect(single_electron.get_prostr_coor(), single_electron.get_dir())
    new_dir = geom.get_new_dir_after_reflect(single_electron.get_dir())
    single_electron.set_coor(new_coor)
    single_electron.set_dir(new_dir)

def transport_process(single_electron, tau_mass, E_mass, dt):

    _make_new_coor(single_electron, dt)
    _make_scatterings(single_electron, tau_mass, E_mass)
    

