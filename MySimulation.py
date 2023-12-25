import numpy as np
import pandas as pd
import Semiconductor as semi

class MySimulation:

    def __init__(self, E_photon, N_eletrons, delta_of_initial_energyes):

        self.E_photon = E_photon
        self.N_eletrons = N_eletrons
        self.delta_of_initial_energyes = delta_of_initial_energyes
