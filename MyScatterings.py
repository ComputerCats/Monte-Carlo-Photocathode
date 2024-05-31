class Semiconductor:

    def __init__(self, E_a, E_g, effective_mass):

        self.E_a = E_a
        self.E_g = E_g
        self.effective_mass = effective_mass

    def get_E_a(self):

        return self.E_a

    def get_E_g(self):

        return self.E_g

    def get_effective_mass(self):

        return self.effective_mass
