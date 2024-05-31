import numpy as np

#Constants

C_CONST = 2.99792458
EV_CONST = 1.602176634
M_E = 9.109

class Electrons:

    def __init__(self, x, y, z, vx, vy, vz):
        
        self.coor = np.array([x, y, z, vx, vy, vz])
        
    def set_electron_propities(self, effective_mass):

        self.effective_mass = effective_mass

    def get_log(self):

        log = f'M_e={self.effective_mass}\nE_a={self.E_a}\nE_g={self.E_g}\ndelta_E_DOS={self.delta_E_DOS}\n' 
    
        return log
    
    def get_E(self):

        return self.effective_mass*1e6*M_E*self.get_veloicity()**2/(2*EV_CONST)

    def set_veloicity(self, new_veloicity):

        self.coor[3:] = new_veloicity

    def add_energy(self, E):

        el_energy = self.get_E()
        
        if el_energy == 0 or el_energy + E < 0:

            raise ValueError(f'Electron energy error, Eel = {el_energy}')

        self.coor[3:] = np.sqrt((el_energy + E)/el_energy)*self.coor[3:]

    def get_veloicity(self):
        
        return np.sqrt(np.dot(self.coor[3:], self.coor[3:]))

    def get_veloicity_vector(self):

        return self.coor[3:]

    def get_effective_mass(self):

        return self.effective_mass

    def get_coor(self):

        return self.coor[:3]

    def set_coor(self, new_coor):

        self.coor[:3] = new_coor

    def add_coor(self, adding_part):

        self.coor[:3] += adding_part




