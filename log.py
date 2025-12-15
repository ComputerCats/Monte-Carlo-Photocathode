import numpy as np
import pandas as pd

class MyLog:

    def __init__(self, name_exp):
        
        self.mass_log = {}
        self.name_exp = name_exp

    def _add_str_to_log(self, param_name, param_value):

        self.mass_log[f'{param_name}'] = param_value

    def save_log(self):

        with open(f'log_{self.name_exp}.txt', 'w') as f:

            for key in self.mass_log:

                f.write(f'{key} = {self.mass_log[key]}\n')

    
