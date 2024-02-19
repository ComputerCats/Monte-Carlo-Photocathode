import numpy as np
import Geometry
import ElTransport as eltrans



def _p_exit(E, E_a, cos_angle):

    E_exit = E*cos_angle*cos_angle

    if E <= 0:

        return 0

    if (np.sqrt(E_a/E) < cos_angle) or (E_exit < E_a):

        return 0

    else:

        result = 4*np.sqrt(E_exit*(E_exit-E_a))/(np.sqrt(E_exit-E_a)+np.sqrt(E_exit))**2

        return result


    return np.random.choice([False, True], p = [1-prop_exit, prop_exit])


    electron_status = geom.get_status(single_electron.get_prostr_coor())

        
            

        else:

            eltrans.reflcation_process(geom, single_electron)


    else:

        return False



            

        



