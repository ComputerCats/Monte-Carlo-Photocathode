import numpy as np

#Constants

C_CONST = 2.99792458
EV_CONST = 1.602176634
M_E = 9.109

def _check_sizes(N_el, mass2):

    if N_el != mass2.shape[0]:

        raise ValueError('Not equal sizes')

class Electrons:

    # [x, y, z, vx, vy, vz, {0 or 1}]
    # last column - "Is electron in simulation (0), if 1: electron exited or died"

    def __init__(self, N_el_in_array: int):

        if N_el_in_array <= 0:

            raise ValueError('N_el_in_array must be greater than zero')

        if N_el_in_array <= 0:

            raise ValueError('N_el_in_array must be grater than zero')

        self.N_el_in_array = N_el_in_array
        
        self.coor = np.zeros((N_el_in_array, 7))

        self.coor[:, -1] = 1

        self._velocity_cache = None
        self._energy_cache = None
        self._cache_valid = False

    def set_velocity(self, new_veloicity, indx_el = None):

        self._cache_valid = False

        if type(indx_el) == type(None):
            _check_sizes(self.N_el_in_array, new_veloicity)
            self.coor[:, 3:-1] = new_veloicity
            return
        if type(indx_el) == int:
            self.coor[indx_el, 3:-1] = new_veloicity
        else:
            self.coor[:, 3:-1][indx_el] = new_veloicity

    def set_electron_properties(self, effective_mass):

        if effective_mass <= 0:

            raise ValueError('Effective mass must be positive')

        self.effective_mass = effective_mass

    def set_coor(self, new_coor, indx_el = None):

        if type(indx_el) == type(None):
            _check_sizes(self.N_el_in_array, new_coor)
            self.coor[:, :-1] = new_coor
            return
        if type(indx_el) == int:
            self.coor[indx_el, :-1] = new_coor

        else:
            self.coor[:, :-1][indx_el] = new_coor

    def set_xcoor(self, new_xcoor, indx_el = None):

        if type(indx_el) == type(None):

            _check_sizes(self.N_el_in_array, new_xcoor)
            self.coor[:, :3] = new_xcoor
            return
        if type(indx_el) == int:

            self.coor[indx_el, :3] = new_xcoor

        else:

            self.coor[:, :3][indx_el] = new_xcoor

    def set_flags(self, new_flags, indx_el = None):
        if type(indx_el) == type(None):
            _check_sizes(self.N_el_in_array, new_flags)
            self.coor[:, -1] = new_flags
            return
        if type(indx_el) == int:

            self.coor[indx_el, -1] = new_flags

        else:

            self.coor[:, -1][indx_el] = new_flags

    def _update_cache(self):

        module_velocity = self.get_module_velocity()
        self._energy_cache = self.effective_mass * 1e6 * M_E * module_velocity**2 / (2 * EV_CONST)
        self._cache_valid = True

    def get_flags(self, indx_el = None):

        if type(indx_el) == type(None):
            return self.coor[:, -1]
        if type(indx_el) == int:
            return self.coor[indx_el, -1]

        else:
            return self.coor[indx_el][:, -1]

    def get_xcoor(self, indx_el = None):

        if type(indx_el) == type(None):
            return self.coor[:, :3]
        if type(indx_el) == int:
            return self.coor[indx_el, :3]
        else:
            return self.coor[indx_el][:, :3]

    def get_log(self):

        log = f'M_e={self.effective_mass}\nE_a={self.E_a}\nE_g={self.E_g}\ndelta_E_DOS={self.delta_E_DOS}\n' 
    
        return log
    
    def get_E(self, indx_el = None):

        if not self._cache_valid:
            self._update_cache()
        if indx_el is None:
            return self._energy_cache
        return self._energy_cache[indx_el]

    def get_N_el(self):

        return self.N_el_in_array

    def get_module_velocity(self, indx_el = None):

        v = self.get_velocity(indx_el)
        
        return np.linalg.norm(v, axis=-1 if v.ndim > 1 else 0)

    def get_velocity(self, indx_el = None):

        if type(indx_el) == type(None):

            return self.coor[:, 3:-1]

        if type(indx_el) == int:
            return self.coor[indx_el, 3:-1]
        else:
            return self.coor[indx_el][:, 3:-1]

    def get_effective_mass(self):

        return self.effective_mass

    def get_coor(self, indx_el = None):

        if type(indx_el) == type(None):
            return self.coor[:, :-1]

        if type(indx_el) == int:
            return self.coor[indx_el, :-1]
        else:
            return self.coor[indx_el][:, :-1]

    def get_alive(self):

        return self.get_flags() == 1

    def add_coor(self, shift, indx_el = None):

        if type(indx_el) == type(None): 
            _check_sizes(self.N_el_in_array, shift)
            self.set_xcoor(self.coor[:, :3] + shift)
            return
        else:
            self.set_xcoor(self.coor[indx_el, :3] + shift, indx_el)

    def add_energy(self, E, indx_el = None):

        el_energy = self.get_E(indx_el)
        new_vel = (np.sqrt((el_energy + E)/el_energy)).reshape(-1, 1)*self.get_velocity(indx_el)
        self.set_velocity(new_vel, indx_el)

    def kill_low_energy_electron(self, kill_energy):
        res = np.where(self.get_E() > kill_energy, self.get_flags(), 0)
        self.set_flags(res)

    def is_end(self):
        return np.all(self.coor[:, -1] == 0)

    def kill_electron(self, indx_el):
        if type(indx_el) == int:
            self.coor[indx_el, -1] = 0
        else:
            self.coor[:, -1][indx_el] = 0

    def is_alive(self, indx_el):
        return self.coor[indx_el, -1] == 1




