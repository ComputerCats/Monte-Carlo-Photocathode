

class Electrons:

    def __init__(self):
        pass

    def set_electron_propities(self, effective_mass):
        self.effective_mass = effective_mass

    def set_semi_params(self, E_g, E_a, delta_E_DOS):
        self.E_a = E_a
        self.E_g = E_g
        self.delta_E_DOS = delta_E_DOS

    def set_transport_params(self, delta_E, tau):
        self.tau = tau
        self.delta_E = delta_E

    def get_log(self):

        log = f'M_e={self.effective_mass}\nE_a={self.E_a}\nE_g={self.E_g}\ndelta_E_DOS={self.delta_E_DOS}\nTau={self.tau}\nDelta_E={self.delta_E}\n' 
    
        return log
    






