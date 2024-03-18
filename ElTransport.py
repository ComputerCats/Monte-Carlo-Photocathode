import numpy as np

# life of electron: (x, y, z, psi, theta, E) -> phonon scattering (change angle, change coor, change energy) -> new iter

def _make_new_coor(single_electron, dt): 

    curr_dir = single_electron.get_dir()
    l_move = single_electron.get_veloicity()*dt

    cos_theta_mass = np.cos(curr_dir[1])
    sin_theta_mass = np.sin(curr_dir[1])
    sin_psi_mass = np.sin(curr_dir[0])
    cos_psi_mass = np.cos(curr_dir[0])

    dx = l_move*sin_theta_mass*cos_psi_mass
    dy = l_move*sin_theta_mass*sin_psi_mass
    dz = l_move*cos_theta_mass

    single_electron.add_coor(np.array([dx, dy, dz]))

def _make_p_l_e_e(single_electron, l_e_e, dt):

    p = 1 - np.exp(-single_electron.get_veloicity()*dt/l_e_e(single_electron.get_E()))

    if p > 1:
        
        raise ValueError('p scat must be less then 1')

    if p < 0:
        
        raise ValueError('p scat must be bigger then 0')

    return p

def _is_scat(p):

    return p > np.random.rand()

def _make_scatterings(single_electron, dt, scatterings_l_e_e, scatterings_E_l_e_e):
    
    new_dir = False

    for indx, l_e in enumerate(scatterings_l_e_e):

        p_scat = _make_p_l_e_e(single_electron, l_e, dt)

        if _is_scat(p_scat):

            single_electron.add_energy(scatterings_E_l_e_e[indx])

            new_dir = True

    if new_dir:

        _make_new_dir(single_electron)

def _make_new_dir(single_electron):

    new_psi = 2*np.pi*np.random.rand()
    new_theta = np.pi*np.random.rand()

    single_electron.set_dir(np.array([new_psi, new_theta]))

def make_initial_dir():

    new_psi = 2*np.pi*np.random.rand()
    new_theta = np.pi*np.random.rand()

    return np.array([new_psi, new_theta])

def transport_process(single_electron, dt, scatterings_l_e_e, scatterings_E_l_e_e):

    _make_new_coor(single_electron, dt)
    _make_scatterings(single_electron, dt, scatterings_l_e_e, scatterings_E_l_e_e)

