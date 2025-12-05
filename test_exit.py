import unittest
import Geometry
import numpy as np
import electron
import ElectronExit as elex
import MyScatterings as scat
import ElTransport as eltrans

#STATUS = {'Exit': 'Exit', 'Outside': 'Outside','Inside': 'Inside','Reflect': 'Reflect', 'Died': 'Died'}
#self.normales = [np.array([1, 0, 0]), np.array([-1, 0, 0]), np.array([0, 0, 1]), np.array([0, 0, -1])]

def get_veloicity_from_energy(E, x, y, z):

    module_vel = 0.001*np.sqrt(E*3.2/(0.11*9.1))

    vec_dir = np.array([x, y, z])

    vec_dir = vec_dir/np.sqrt(np.dot(vec_dir, vec_dir))

    result = module_vel*vec_dir

    return result

class TestConvexShape(unittest.TestCase):

    def setUp(self):

        plane_cathode = plane_cathode = Geometry.Plane(np.array([0, 0, 0.512]), np.array([0, 0, -1]))
        self.geom = Geometry.ConvexShape()
        self.geom.add_plane(plane_cathode, True)
        self.semiconductor = scat.Semiconductor(0.6, 1.1, 0.11)

    def test_exit_process(self):

        kill_energy = 0.3

        single_electron1 = electron.Electrons(0, 0, 0.51, 0, 0, -1)
        single_electron1.set_electron_propities(0.11)
        vel = get_veloicity_from_energy(1.8, 0, 0, -1)

        single_electron1.set_veloicity(vel)

        self.assertEqual(elex.exit_process(self.geom, single_electron1, self.semiconductor, kill_energy), elex.EXIT_PROCESS_STATUS['Out'])

        ###

        single_electron1 = electron.Electrons(0, 0, 0.51, 0, 0, 1)
        single_electron1.set_electron_propities(0.11)
        vel = get_veloicity_from_energy(1.8, 0, 0, 1)

        single_electron1.set_veloicity(vel)

        self.assertEqual(elex.exit_process(self.geom, single_electron1, self.semiconductor, kill_energy), elex.EXIT_PROCESS_STATUS['New_iter'])

        ###

        single_electron1 = electron.Electrons(0, 0, 0.51, 0, 1, -1)
        single_electron1.set_electron_propities(0.11)
        vel = get_veloicity_from_energy(1.8, 0, 1, 0)

        single_electron1.set_veloicity(vel)

        self.assertEqual(elex.exit_process(self.geom, single_electron1, self.semiconductor, kill_energy), elex.EXIT_PROCESS_STATUS['New_iter'])

        ###

        single_electron1 = electron.Electrons(0, 0, 0.51, 0, 0, -1)
        single_electron1.set_electron_propities(0.11)
        vel = get_veloicity_from_energy(0.4, 0, 0, -1)

        single_electron1.set_veloicity(vel)

        self.assertEqual(elex.exit_process(self.geom, single_electron1, self.semiconductor, kill_energy), elex.EXIT_PROCESS_STATUS['New_iter'])

        ###

        single_electron1 = electron.Electrons(0, 0, 0.51, 0, 0, -1)
        single_electron1.set_electron_propities(0.11)
        vel = get_veloicity_from_energy(1.8, 1, 0, 0)

        single_electron1.set_veloicity(vel)

        self.assertEqual(elex.exit_process(self.geom, single_electron1, self.semiconductor, kill_energy), elex.EXIT_PROCESS_STATUS['New_iter'])

        ###

        single_electron1 = electron.Electrons(0, 0, 0.51, 0, 0, -1)
        single_electron1.set_electron_propities(0.11)
        vel = get_veloicity_from_energy(2, 0, 0.1, -0.9)
        single_electron1.set_veloicity(vel)

        self.assertEqual(elex.exit_process(self.geom, single_electron1, self.semiconductor, kill_energy), elex.EXIT_PROCESS_STATUS['Out'])

    def test_reflectrion(self):

        kill_energy = 0.3

        single_electron1 = electron.Electrons(0, 0, 0.51, 0, 0, -1)
        single_electron1.set_electron_propities(0.11)
        vel = get_veloicity_from_energy(0.8, 0, 0.5, -0.5)

        single_electron1.set_veloicity(np.copy(vel))
        elex.exit_process(self.geom, single_electron1, self.semiconductor, kill_energy)

        vel[2] = -vel[2]
        print(single_electron1.get_veloicity())
        print(vel)
        self.assertEqual(single_electron1.get_veloicity() == vel, True)
        

if __name__ == '__main__':
    unittest.main()