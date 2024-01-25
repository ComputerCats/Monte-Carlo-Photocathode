import numpy as np

class Electrons:

    def __init__(self, x, y, z, psi, theta, E):
        
        self.coor = np.array([x, y, z, psi, theta, E])

    def set_electron_propities(self, effective_mass):

        self.effective_mass = effective_mass

    def get_log(self):

        log = f'M_e={self.effective_mass}\nE_a={self.E_a}\nE_g={self.E_g}\ndelta_E_DOS={self.delta_E_DOS}\n' 
    
        return log
    
    def get_E(self):

        return self.coor[-1]

    def add_energy(self, E):

        self.coor[-1] += E

    def set_E(self, new_energy):

        self.coor[-1] = new_energy

    def get_veloicity(self):

        return 0.001*np.sqrt(self.coor[-1]*2*1.6/(self.effective_mass*9.1))

    def get_prostr_coor(self):

        return self.coor[:3]

    def set_coor(self, new_coor):

        self.coor[:3] = new_coor

    def add_coor(self, adding_part):

        self.coor[:3] += adding_part

    def get_dir(self):

        return self.coor[3:5]

    def set_dir(self, new_dir):

        self.coor[3:5] = new_dir


