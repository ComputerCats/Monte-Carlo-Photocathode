import numpy as np
import Geometry
import ElTransport as eltrans

#Geometry.STATUS = {'Exit': 'Exit', 'Outside': 'Outside','Inside': 'Inside','Reflect': 'Reflect', 'Died': 'Died'}

EXIT_PROCESS_STATUS = {'Out': 'Out', 'Died': 'Died', 'New_iter': 'New_iter'}

def _reflcation_process(geom, single_electron):

    new_coors = geom.get_new_coors_after_reflect(single_electron)
    single_electron.set_coor(new_coors['new_coor'])
    single_electron.set_veloicity(new_coors['new_dir'])

def _p_exit(E, E_a, cos_angle):
    
    E_exit = E*cos_angle*cos_angle
    
    if E <= 0:

        return 0

    if (np.sqrt(E_a/E) < cos_angle) and (E_exit > E_a):

        result = 4*np.sqrt(E_exit*(E_exit-E_a))/(np.sqrt(E_exit-E_a)+np.sqrt(E_exit))**2

        return result

    else:

        return 0

def _is_exit(geom, single_electron, semiconductor):

    prop_exit = _p_exit(single_electron.get_E(), semiconductor.get_E_a(), geom.get_cos_angle(single_electron))
    
    return prop_exit > np.random.rand()

def _is_low_energy_electron(single_electron, kill_energy):

        if single_electron.get_E() < kill_energy:

            return True

        else:

            return False

def exit_process(geom, single_electron, semiconductor, kill_energy):
    
    electron_status = geom.get_status(single_electron)

    if _is_low_energy_electron(single_electron, kill_energy):   # low energy case
        
        return EXIT_PROCESS_STATUS['Died']

    if electron_status == Geometry.STATUS['Exit']:              # exit process
        
        if _is_exit(geom, single_electron, semiconductor):
            
            return EXIT_PROCESS_STATUS['Out']

        else:

            _reflcation_process(geom, single_electron)

            return EXIT_PROCESS_STATUS['New_iter']

    if electron_status == Geometry.STATUS['Died']:              

        return EXIT_PROCESS_STATUS['Died']

    if electron_status == Geometry.STATUS['Inside']:
        
        return EXIT_PROCESS_STATUS['New_iter']

    if electron_status == Geometry.STATUS['Reflect']:
        print('Died')
        _reflcation_process(geom, single_electron)

        return EXIT_PROCESS_STATUS['New_iter']


            

        



