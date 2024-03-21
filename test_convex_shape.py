import unittest
import Geometry
import numpy as np
import electron

#STATUS = {'Exit': 'Exit', 'Outside': 'Outside','Inside': 'Inside','Reflect': 'Reflect', 'Died': 'Died'}
#self.normales = [np.array([1, 0, 0]), np.array([-1, 0, 0]), np.array([0, 0, 1]), np.array([0, 0, -1])]

class TestConvexShape(unittest.TestCase):

    def setUp(self):

        plane_cathode = plane_cathode = Geometry.Plane(np.array([0, 0, 0.512]), np.array([0, 0, -1]))
        self.geom = Geometry.ConvexShape()
        self.geom.add_plane(plane_cathode, True)

    def test_status(self):

        single_electron1 = electron.Electrons(0, 0, 0.510, 0, 0, 0)
        self.assertEqual(self.geom.get_status(single_electron1), Geometry.STATUS['Exit'])

        single_electron2 = electron.Electrons(0, 0, 0.514, 0, 0, 0)
        self.assertEqual(self.geom.get_status(single_electron2), Geometry.STATUS['Inside'])

        single_electron3 = electron.Electrons(0, 2, 0.6, 0, 0, 0)
        self.assertEqual(self.geom.get_status(single_electron3), Geometry.STATUS['Inside'])

        single_electron4 = electron.Electrons(0, -2, -0.6, 0, 0, 0)
        self.assertEqual(self.geom.get_status(single_electron4), Geometry.STATUS['Exit'])

        single_electron5 = electron.Electrons(0, -1, 0.512, 0, 0, 0)
        self.assertEqual(self.geom.get_status(single_electron5), Geometry.STATUS['Inside'])
    
    def test_get_cos_angle(self):

        single_electron1 = electron.Electrons(0, 0, 0.51, 0, 0, -1)
        self.assertEqual(self.geom.get_cos_angle(single_electron1), 1)

        single_electron2 = electron.Electrons(0, 0, 0.51, 0, 0, 1)
        self.assertEqual(self.geom.get_cos_angle(single_electron2), -1)

        single_electron3 = electron.Electrons(0, 2, 0.51, 0, 1, 0)
        self.assertEqual(self.geom.get_cos_angle(single_electron3), 0)

        single_electron4 = electron.Electrons(0, -2, 0.51, 0, -1, 0)
        self.assertEqual(self.geom.get_cos_angle(single_electron4), 0)

        single_electron5 = electron.Electrons(0, -1, 0.51, 0, 1, -1)
        self.assertEqual(round(self.geom.get_cos_angle(single_electron5), 3), round(np.cos(np.pi/4), 3))

        single_electron6 = electron.Electrons(0, -1, 0.51, 0, 1, 1)
        self.assertEqual(round(self.geom.get_cos_angle(single_electron6), 3), -round(np.cos(np.pi/4), 3))

        single_electron7 = electron.Electrons(0, -1, 0.51, 1, 0, -1)
        self.assertEqual(round(self.geom.get_cos_angle(single_electron7), 3), round(np.cos(np.pi/4), 3))
    

if __name__ == '__main__':
    unittest.main()
