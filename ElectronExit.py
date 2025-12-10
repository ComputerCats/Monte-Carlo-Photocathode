import numpy as np
import Geometry
import ElTransport as eltrans

def _p_exit(E, E_a, cos_angle):
    
    E_exit = E*cos_angle*cos_angle

    cond_out = (np.sqrt(E_a/E) < cos_angle) * (E_exit > E_a)
    
    result = cond_out * 4*np.sqrt(E_exit*(np.abs(E_exit-E_a)))/(np.sqrt(np.abs(E_exit-E_a))+np.sqrt(E_exit))**2
    
    return result

def is_exit(geom, electrons, semiconductor, indx_el: int, indx_out: int) -> bool:

    prop_exit = _p_exit(electrons.get_E(indx_el), semiconductor.get_E_a(), geom.get_cos_angle(electrons, indx_el, indx_out))
    
    return prop_exit > np.random.rand()

def _get_this_electron_emitance(geom, electrons, el_indx: int, out_plane_indx: int):

    cos_out = geom.get_cos_angle(electrons, el_indx, out_plane_indx)

    return electrons.get_E(el_indx)*(1-cos_out*cos_out)

def _reflection_process(geom, electrons, indx_el, indx_plane):

    new_coor = geom.get_new_coors_after_reflect(electrons, indx_el, indx_plane)

    electrons.set_coor(new_coor, indx_el)

def exit_process(geom, electrons, semiconductor):

    N_exit  = 0
    emmitance = 0

    N_el_sub_sim = electrons.get_N_el()

    for indx_el in range(N_el_sub_sim):

        if electrons.is_alive(indx_el):
            indx_out_plane = geom.get_out_plane(electrons, indx_el)

            if indx_out_plane > -1: # refactor it!

                type_exit_plane = geom.get_type_exit_plane(indx_out_plane)
                if type_exit_plane:
                    if is_exit(geom, electrons, semiconductor, indx_el, indx_out_plane):
                        N_exit += 1
                        emmitance += _get_this_electron_emitance(geom, electrons, indx_el, indx_out_plane)
                        electrons.kill_electron(indx_el)
                    else:
                        _reflection_process(geom, electrons, indx_el, indx_out_plane)

                else:
                    electrons.kill_electron(indx_el)

    return {'N_exit': N_exit, 'Emmitance': emmitance}

            

        



