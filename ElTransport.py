import numpy as np

# life of electron: (x, y, z, psi, theta, E) -> phonon scattering (change angle, change coor, change energy) -> new iter

C_CONST = 2.99792458
EV_CONST = 1.602176634
M_E = 9.109

def _make_new_coor(electrons, dt): 

    l_move = electrons.get_velocity()*dt

    electrons.add_coor(l_move)

def _make_p_l_e_e(electrons, indx_el, l_e_e, dt):

    p = 1 - np.exp(-electrons.get_module_velocity(indx_el)*dt/l_e_e(electrons.get_E(indx_el)))

    return p

def _is_scat(p):

    return p > np.random.rand()

def _make_scatterings(electrons, indx_el, dt, scatterings_l_e_e, scatterings_E_l_e_e):
    
    new_dir = False

    for indx, l_e in enumerate(scatterings_l_e_e):

        p_scat = _make_p_l_e_e(electrons, indx_el, l_e, dt)

        if _is_scat(p_scat):

            if electrons.get_E(indx_el) + scatterings_E_l_e_e[indx] <= 0.01: break #????

            electrons.add_energy(scatterings_E_l_e_e[indx], indx_el)

            new_dir = True

    if new_dir:

        _make_new_dir(electrons, indx_el)

def _make_new_dir(electrons, el_indx):

    new_psi = 2*np.pi*np.random.rand()
    new_theta = np.pi*np.random.rand()
    
    module_vel = electrons.get_module_velocity(el_indx)

    vx = module_vel*np.sin(new_theta)*np.cos(new_psi)
    vy = module_vel*np.sin(new_theta)*np.sin(new_psi)
    vz = module_vel*np.cos(new_theta)

    electrons.set_velocity(np.array([vx, vy, vz]), el_indx)

def make_initial_veloicity(electrons, energy):

    vec = np.random.normal(0, 1, (electrons.get_N_el_in_ar(), 3))
    
    # Нормализуем каждый вектор к длине 1
    norms = np.linalg.norm(vec, axis=1, keepdims=True)
    unit_vectors = vec / norms

    module_vel = 1e-3*np.sqrt(energy*2*EV_CONST/(electrons.get_effective_mass()*M_E))

    electrons.set_velocity(module_vel*unit_vectors)

def transport_process(electrons, dt, scatterings_l_e_e, scatterings_E_l_e_e):

    _make_new_coor(electrons, dt)
    
    for indx_el in range(electrons.get_N_el_in_ar()):
        _make_scatterings(electrons, indx_el, dt, scatterings_l_e_e, scatterings_E_l_e_e)

