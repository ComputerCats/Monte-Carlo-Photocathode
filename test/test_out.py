import ElectronExit as exit
import electron as el
import numpy as np
import Geometry as geom
import MyScatterings as scat

def test_exit_1():

    halfscpace = geom.HalfspaceGeom(np.array([0, 0, 0]), np.array([0, 0, 1]))

    single_electron = el.Electrons(0, 0, -0.05, 0, np.pi, 1)
    print(np.cos(single_electron.get_dir()[1]))

    semiconductor = scat.Semiconductor(0.3, 1.6, 0.12)
    print('out test')
    print(f'test_1: {exit.exit_process(halfscpace, single_electron, semiconductor)}, pass = {exit.exit_process(halfscpace, single_electron, semiconductor) == True}')

def test_exit_2():

    halfscpace = geom.HalfspaceGeom(np.array([0, 0, 0]), np.array([0, 0, 1]))

    single_electron = el.Electrons(0, 0, 0.05, 0, np.pi, 1)
    print(np.cos(single_electron.get_dir()[1]))

    semiconductor = scat.Semiconductor(0.3, 1.6, 0.12)

    print('inside test')
    print(f'test_2: {exit.exit_process(halfscpace, single_electron, semiconductor)}, pass = {exit.exit_process(halfscpace, single_electron, semiconductor) == False}')

def test_exit_3():

    halfscpace = geom.HalfspaceGeom(np.array([0, 0, 0]), np.array([0, 0, 1]))

    single_electron = el.Electrons(0, 0, -0.05, 0, np.pi-np.pi/3-np.pi/12, 1)
    print(np.cos(single_electron.get_dir()[1]))

    semiconductor = scat.Semiconductor(0.3, 1.6, 0.12)
    print('reflection test')
    print(f'test_3: {exit.exit_process(halfscpace, single_electron, semiconductor)}, pass = {exit.exit_process(halfscpace, single_electron, semiconductor) == False}')

def test_exit_4():

    halfscpace = geom.HalfspaceGeom(np.array([0, 0, 0]), np.array([0, 0, 1]))

    single_electron = el.Electrons(0, 0, -0.05, 0, np.pi-np.pi/3-np.pi/12, 1)
    print(np.cos(single_electron.get_dir()[1]))

    semiconductor = scat.Semiconductor(0.3, 1.6, 0.12)
    exit.exit_process(halfscpace, single_electron, semiconductor)
    print('reflection movement test')
    print(f'test_4: {single_electron.get_dir()}, pass = {single_electron.get_dir()[1]} == {np.pi/3+np.pi/12}')
    print(f'test_4: {single_electron.get_prostr_coor()} pass = {single_electron.get_prostr_coor()[2] == 0.05}')

def test_exit_5():

    halfscpace = geom.HalfspaceGeom(np.array([0, 0, 0]), np.array([0, 0, 1]))

    single_electron = el.Electrons(0, 0, -0.05, 0, np.pi-np.pi/3-np.pi/12, 0.2)
    print(np.cos(single_electron.get_dir()[1]))

    semiconductor = scat.Semiconductor(0.3, 1.6, 0.12)
    result = exit.exit_process(halfscpace, single_electron, semiconductor)
    print('low energy test')
    print(f'test_5: {result}, pass = {result == False}')


test_exit_1()
test_exit_2()
test_exit_3()
test_exit_4()
test_exit_5()