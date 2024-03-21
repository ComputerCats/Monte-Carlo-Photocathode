import numpy as np

# life of electron: (x, y, z, psi, theta, E) -> phonon scattering (change angle, change coor, change energy) -> new iter

def _make_new_coor(single_electron, dt): 

    l_move = single_electron.get_veloicity_vector()*dt

    single_electron.add_coor(l_move)

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

    module_vel = single_electron.get_veloicity()

    vx = module_vel*np.sin(new_theta)*np.cos(new_psi)
    vy = module_vel*np.sin(new_theta)*np.sin(new_psi)
    vz = module_vel*np.cos(new_theta)

    single_electron.set_veloicity(np.array([vx, vy, vz]))

def make_initial_dir(single_electron, energy):

    new_psi = 2*np.pi*np.random.rand()
    new_theta = np.pi*np.random.rand()

    module_vel = 0.001*np.sqrt(energy*3.2/(single_electron.get_effective_mass()*9.1))

    vx = module_vel*np.sin(new_theta)*np.cos(new_psi)
    vy = module_vel*np.sin(new_theta)*np.sin(new_psi)
    vz = module_vel*np.cos(new_theta)

    new_veloicity = np.array([vx, vy, vz])

    return new_veloicity

def transport_process(single_electron, dt, scatterings_l_e_e, scatterings_E_l_e_e):

    _make_new_coor(single_electron, dt)
    _make_scatterings(single_electron, dt, scatterings_l_e_e, scatterings_E_l_e_e)

