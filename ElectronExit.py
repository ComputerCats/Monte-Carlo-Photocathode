import numpy as np
import Geometry
import ElTransport as eltrans

#Geometry.STATUS = {'Exit': 'Exit','Inside': 'Inside','Reflect': 'Reflect','~reflect': '~reflect'}

def p_exit(E, E_a, cos_angle):

    E_exit = E*cos_angle*cos_angle

    if E <= 0:

        return 0

    if (np.sqrt(E_a/E) < cos_angle) or (E_exit < E_a):

        return 0

    else:

        result = 4*np.sqrt(E_exit*(E_exit-E_a))/(np.sqrt(E_exit-E_a)+np.sqrt(E_exit))**2

        return result

def is_exit(geom, single_electron, semiconductor):

    prop_exit = p_exit(single_electron.get_E(), semiconductor.get_E_a(), np.cos(single_electron.get_dir()[1]))
    return np.random.choice([False, True], p = [1-prop_exit, prop_exit])

def exit_process(geom, single_electron, semiconductor):

    electron_status = geom.get_status(single_electron)

    if electron_status == Geometry.STATUS['Exit']:
        
        if is_exit(geom, single_electron, semiconductor):
            
            return True

        else:

            eltrans.reflcation_process(geom, single_electron)

            return False

    else:

        return False



            

        



