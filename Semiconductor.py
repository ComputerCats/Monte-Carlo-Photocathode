import numpy as np
import pandas as pd
import MySimulation as MySim
import electron as electr


class Semiconductor:

    def __init__(self, E_g, E_a, delta_phonon):

        self.E_g = E_g
        self.E_a = E_a
        self.delta_phonon = delta_phonon