import numpy as np
import Geometry

#Geometry.STATUS = {'Exit': 'Exit','Inside': 'Inside','Reflect': 'Reflect','~reflect': '~reflect'}

def p_exit(E, E_a, cos_angle):

    E_exit = E*cos_angle*cos_angle

    if E <= 0:

        return 0

    if (np.sqrt(E_a/E) < cos_angle) and (E_exit > E_a):

        result = 4*np.sqrt(E_exit*(E_exit-E_a))/(np.sqrt(E_exit-E_a)+np.sqrt(E_exit))**2

        return result

    else:

        return 0

def _is_exit(geom, single_electron):

    prop_exit = p_exit(single_electron.get_E(), single_electron.E_a, geom.get_cos_angle(single_electron))

    return np.random.choice([False, True], p = [1-prop_exit, prop_exit])

def exit_process(geom, single_electron):

    status = geom.get_status(single_electron)

    if status == Geometry.STATUS['Exit']:

        if _is_exit(geom, single_electron):

            return True

        else:

            geom.reflect(single_electron)

            return False

    else:

        return False



