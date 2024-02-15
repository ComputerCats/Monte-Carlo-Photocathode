import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

sys.path.append(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')

import electron

def test_make_new_dir(way_to):

    #test dir distribution
    initial_N_electrons = 10
    electron_gas = np.zeros((initial_N_electrons, 5))
    dirs = electron.make_new_dir(electron_gas)

    print('dirs = ' + str(dirs[:, 3:5]))

test_make_new_dir(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')