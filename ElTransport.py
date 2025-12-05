import numpy as np

# life of electron: (x, y, z, vx, vy, vz, {0, 1}) -> phonon scattering (change velocity (if needed), change coor) -> new iter

C_CONST = 2.99792458
EV_CONST = 1.602176634
M_E = 9.109

def _make_new_coor(electrons, dt): 

    indx_alive = electrons.get_flags() < 1

    l_move = electrons.get_velocity(indx_alive)*dt

    electrons.add_coor(l_move, indx_alive)

def _make_p_l_e_e(electrons, l_e_e, dt):

    p = 1 - np.exp(-electrons.get_module_velocity()*dt/l_e_e(electrons.get_E()))

    return p

def _is_scat(p):

    return p > np.random.rand()

def _make_scatterings(electrons, dt, scatterings_l_e_e, scatterings_E): # TEST IT
    
    n_electrons = electrons.get_N_el()
    n_scat = len(scatterings_l_e_e)

    mass_el_scat = np.zeros((n_electrons, n_scat)).astype(bool)

    for indx, l_e in enumerate(scatterings_l_e_e):

        mass_el_scat[:, indx] = _make_p_l_e_e(electrons, l_e, dt) > np.random.rand(1, n_electrons)

    any_scattering = ((mass_el_scat*(electrons.get_flags()).astype(bool).reshape(-1, 1)).any(axis=1))
    
    if np.any(any_scattering):

        N_el_have_scat = len(np.where(any_scattering)[0])
        indx_scat = np.where(any_scattering)[0]
        new_energy = np.zeros((N_el_have_scat, 1))

        for indx, scat_en in enumerate(scatterings_E):
            new_energy[:, 0] += scat_en*mass_el_scat[any_scattering][:, indx]

        _make_velocities_after_scattering(electrons, new_energy, any_scattering, N_el_have_scat)

def _make_velocities_after_scattering(electrons, delta_E, any_scattering, N_el_have_scat):

    vec = np.random.normal(0, 1, (N_el_have_scat, 3))
    
    norms = np.linalg.norm(vec, axis=1, keepdims=True)
    unit_vectors = vec / norms
    energy = electrons.get_E(any_scattering).reshape(-1, 1) + delta_E
    module_vel = 1e-3*np.sqrt(energy*2*EV_CONST/(electrons.get_effective_mass()*M_E))
    electrons.set_velocity(module_vel*unit_vectors, any_scattering)

def make_initial_veloicity(electrons, energy):

    vec = np.random.normal(0, 1, (electrons.get_N_el(), 3))
    
    norms = np.linalg.norm(vec, axis=1, keepdims=True)
    unit_vectors = vec / norms

    module_vel = 1e-3*np.sqrt(energy*2*EV_CONST/(electrons.get_effective_mass()*M_E))

    electrons.set_velocity(module_vel*unit_vectors)

def transport_process(electrons, dt, scatterings_l_e_e, scatterings_E):

    _make_new_coor(electrons, dt)
    _make_scatterings(electrons, dt, scatterings_l_e_e, scatterings_E)

