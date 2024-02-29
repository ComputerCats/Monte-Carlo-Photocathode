import multiprocessing as mp
import numpy as np
import MonteCarlo as monte

class ParallelMonteCarlo:

    def set_parallel_params(self, N_kers, N_electrons):

        self.N_kers = N_kers
        self.N_electrons = N_electrons

    def set_params(self, sim_params):

        self.sim_params = sim_params

    def _init_params(self):

        self.sims = []

    def run_parallel(self):

        with mp.Pool(self.N_kers) as p:
            super().run_simulation()
