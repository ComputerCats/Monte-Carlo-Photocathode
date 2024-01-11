import numpy as np
import matplotlib.pyplot as plt
import sys
import Distributions
import spectrum_maker_Cs3Sb as spec


sys.path.append(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')


def test_make_coordinate_DOS(way_to):

    #halfspace coor distribution, 1d distribution
    coor_DOS = Distributions.make_coordinate_DOS(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\fdtd\halfspace\V_detector_f0010.d')
    print(f'result_dos = {coor_DOS}')
    print(f'norm must be equal 1, curr sum = {np.cumsum(coor_DOS[:, -1], axis = 0)[-1]}')

    plt.plot(coor_DOS[:, 2], coor_DOS[:, 3])
    plt.grid()
    plt.xlabel('z')
    plt.ylabel('DOS')
    plt.savefig('coor_dos.png')

def test_make_coordinate_DOS_from_alpha():

    energies = np.linspace(1.95, 2.5, 5)

    fig, ax = plt.subplots()

    for energy in energies:

        coor_dos = spec.get_coor_Cs3Sb(energy)

        ax.plot(coor_dos[:, 2], 100*coor_dos[:, 3], label = f'$\gamma$ = {energy}')

    ax.grid()
    ax.set_xlabel('$\Delta$z, $\mu$')
    ax.set_ylabel('DOS, %')
    ax.legend()
    fig.savefig('distr/Absorption_Cs3Sb_energy.png')

#test_make_coordinate_DOS(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')
test_make_coordinate_DOS_from_alpha()