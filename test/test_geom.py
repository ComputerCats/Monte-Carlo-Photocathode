import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import Geometry

sys.path.append(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')

import electron


def test_geom(way_to):

    way_to = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode'

    geom = Geometry.HalfspaceGeom(np.array([0, 0, 0.512]), np.array([0, 0, 1]))
    point_1 = np.array([0, 0, 0.53]) #inside
    point_2 = np.array([0, 0, 2]) #inside
    point_3 = np.array([0, 0, -0.2]) #Exit
    point_4 = np.array([0, 0, 0.511]) #Exit
    point_5 = np.array([0, 0, 0.513]) #inside

    status1 = geom.get_status(point_1)
    status2 = geom.get_status(point_2)
    status3 = geom.get_status(point_3)
    status4 = geom.get_status(point_4)
    status5 = geom.get_status(point_5)

    if status1 == 'Inside':
        print(status1)
        print('test1: pass')

    else:
        print(status1)
        print('test1: not pass')

    if status2 == 'Inside':
        print(status2)
        print('test2: pass')

    else:
        print(status2)
        print('test2: not pass')

    if status3 == 'Exit':
        print(status3)
        print('test3: pass')

    else:
        print(status3)
        print('test3: not pass')

    if status4 == 'Exit':
        print(status4)
        print('test4: pass')

    else:
        print(status4)
        print('test4: not pass')

    if status5 == 'Inside':
        print(status5)
        print('test5: pass')

    else:
        print(status5)
        print('test5: not pass')


test_geom(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode')