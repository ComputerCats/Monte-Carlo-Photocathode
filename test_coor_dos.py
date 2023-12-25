import numpy as np
import matplotlib.pyplot as plt
import sys
import Distributions

sys.path.append(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')


def test_make_coordinate_DOS(way_to):

    #halfspace coor distribution, 1d distribution
    coor_DOS = Distributions.make_coordinate_DOS(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\fdtd\halfspace\V_detector_f0010.d')
    print(f'result_dos = {coor_DOS}')
    print(f'norm must be equal 1, curr sum = {np.cumsum(coor_DOS[:, -1], axis = 0)[-1]}')

    plt.plot(coor_DOS[:, 2], coor_DOS[:, -1])
    plt.grid()
    plt.xlabel('z')
    plt.ylabel('DOS')
    plt.savefig('coor_dos.png')

test_make_coordinate_DOS(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')
