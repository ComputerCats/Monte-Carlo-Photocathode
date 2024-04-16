import unittest
import Geometry
import numpy as np
import electron

STATUS = {'Exit': 'Exit', 'Outside': 'Outside','Inside': 'Inside','Reflect': 'Reflect', 'Died': 'Died'}

class TestConvexShape(unittest.TestCase):

    def setUp(self):

        ###PLANE

        plane_cathode = Geometry.Plane(np.array([0, 0, 0.512]), np.array([0, 0, -1]))
        self.geom = Geometry.ConvexShape()
        self.geom.add_plane(plane_cathode, STATUS['Exit'])

        self.h = 0.04
        self.z0 = 0.512

        ###CATHODE

        self.geom_cathode = Geometry.ConvexShape()

        plane_cathode = Geometry.Plane(np.array([0, 0, self.z0-self.h]), np.array([0, 0, -1]))
        plane_down = Geometry.Plane(np.array([0, 0, self.z0]), np.array([0, 0, 1]))

        self.geom_cathode.add_plane(plane_cathode, STATUS['Exit'])
        self.geom_cathode.add_plane(plane_down, STATUS['Died'])

        ###BOX

        self.geom_box = Geometry.ConvexShape()

        self.h_cathode = 0.024
        self.w = 0.2

        point_1 = np.array([0, 0, 0])
        point_2 = np.array([self.w, 0, -self.h_cathode])

        normale_box_up_up = np.array([0, 0, -1])
        normale_box_up_right = np.array([1, 0, 0])
        normale_box_down_left = np.array([-1, 0, 0])
        normale_box_down_up = np.array([0, 0, 1])

        plane_cathode_up_up = Geometry.Plane(point_2, normale_box_up_up)
        plane_cathode_up_right = Geometry.Plane(point_2, normale_box_up_right)
        plane_cathode_down_left = Geometry.Plane(point_1, normale_box_down_left)
        plane_cathode_down_up = Geometry.Plane(point_1, normale_box_down_up)

        self.geom_box.add_plane(plane_cathode_up_up, STATUS['Exit'])
        self.geom_box.add_plane(plane_cathode_up_right, STATUS['Exit'])
        self.geom_box.add_plane(plane_cathode_down_left, STATUS['Exit'])
        self.geom_box.add_plane(plane_cathode_down_up, STATUS['Died'])

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

    def test_cathode(self):

        single_electron1 = electron.Electrons(0, 0, self.z0-self.h - 0.002, 0, 0, 0)
        self.assertEqual(self.geom_cathode.get_status(single_electron1), Geometry.STATUS['Exit'])

        single_electron2 = electron.Electrons(self.w/2, 0, self.z0-self.h/2, 0, 0, 0)
        self.assertEqual(self.geom_cathode.get_status(single_electron2), Geometry.STATUS['Inside'])

        single_electron3 = electron.Electrons(self.w/2, 0, 0.002, 0, 0, 0)
        self.assertEqual(self.geom_cathode.get_status(single_electron3), Geometry.STATUS['Exit'])

        single_electron4 = electron.Electrons(-0.002, 0, self.z0+0.002, 0, 0, 0)
        self.assertEqual(self.geom_cathode.get_status(single_electron4), Geometry.STATUS['Died'])

    def test_box(self):

        single_electron1 = electron.Electrons(0, 0, -self.h_cathode - 0.002, 0, 0, 0)
        self.assertEqual(self.geom_box.get_status(single_electron1), Geometry.STATUS['Exit'])

        single_electron2 = electron.Electrons(self.w/2, 0, -2*self.h_cathode/3, 0, 0, 0)
        self.assertEqual(self.geom_box.get_status(single_electron2), Geometry.STATUS['Inside'])

        single_electron3 = electron.Electrons(self.w/2, 2, 0.005, 0, 0, 0)
        self.assertEqual(self.geom_box.get_status(single_electron3), Geometry.STATUS['Died'])

        single_electron4 = electron.Electrons(0, -2, -0.6, 0, 0, 0)
        self.assertEqual(self.geom_box.get_status(single_electron4), Geometry.STATUS['Exit'])

        single_electron5 = electron.Electrons(0, -1, -self.h_cathode/2, 0, 0, 0)
        self.assertEqual(self.geom_box.get_status(single_electron5), Geometry.STATUS['Inside'])

        single_electron6 = electron.Electrons(0, -1, 0.555, 0, 0, 0)
        self.assertEqual(self.geom_box.get_status(single_electron6), Geometry.STATUS['Died'])

        single_electron7 = electron.Electrons(0, -1, 0.6, 0, 0, 0)
        self.assertEqual(self.geom_box.get_status(single_electron7), Geometry.STATUS['Died'])

        single_electron8 = electron.Electrons(0, -1, 0.513, 0, 0, 0)
        self.assertEqual(self.geom_box.get_status(single_electron8), Geometry.STATUS['Died'])
    
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
