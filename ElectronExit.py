import numpy as np
import Geometry
import ElTransport as eltrans

def _p_exit(E, E_a, cos_angle):
    
    E_exit = E*cos_angle*cos_angle

    cond_out = (np.sqrt(E_a/E) < cos_angle) * (E_exit > E_a)

    result = cond_out * 4*np.sqrt(E_exit*(np.abs(E_exit-E_a)))/(np.sqrt(np.abs(E_exit-E_a))+np.sqrt(E_exit))**2

    return result

def is_exit(geom, electrons, semiconductor, indx_el: int) -> int:

    prop_exit = _p_exit(electrons.get_E(indx_el), semiconductor.get_E_a(), geom.get_cos_angle(electrons, indx_el))
    
    return prop_exit > np.random.rand()

def _is_not_low_energy_electron(electrons, kill_energy):

        return electrons.get_E() > kill_energy

def _get_this_electron_emitance(self, electrons, el_indx: int, out_plane_indx: int):

    cos_out = self.geometry.get_cos_angle(electrons, el_indx, out_plane_indx)

    return electrons.get_E(el_indx)*(1-cos_out*cos_out)

def exit_process(geom, electrons, semiconductor, kill_energy, N_exit, emmitance):
    
    # electron_status = 1 if Exit
    # electron_status = 0 if Inside

    electrons.set_electron_flags(electrons.get_electron_flags()*_is_not_low_energy_electron(electrons, kill_energy))   # low energy case

    N_el_sub_sim = electrons.get_N_el_in_ar()

    for indx_el in range(N_el_sub_sim):

        if not electrons.get_flags(indx_el):

            indx_out_plane = geom.get_out_plane()

            if indx_out_plane > -1: # refactor it!

                type_exit_plane = geom.get_type_exit_plane(indx_out_plane)
                if type_exit_plane:
                    if is_exit(geom, electrons, semiconductor, indx_el):

                        N_exit += type_exit_plane*1
                        emmitance += type_exit_plane*_get_this_electron_emitance(electrons, indx_el, indx_out_plane)

                else:
                    electrons.kill_electron(indx_el)

    return {'N_exit': N_exit, 'Emmitance': emmitance}

            

        



