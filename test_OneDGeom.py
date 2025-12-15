import unittest
import Geometry
import numpy as np

#STATUS = {'Exit': 'Exit', 'Outside': 'Outside','Inside': 'Inside','Reflect': 'Reflect', 'Died': 'Died'}
#self.normales = [np.array([1, 0, 0]), np.array([-1, 0, 0]), np.array([0, 0, 1]), np.array([0, 0, -1])]

class TestOneDGeom(unittest.TestCase):

    def setUp(self):
        
        self.pos_0 = np.array([0, 0, 0.512])
        self.h = 0.032
        self.h_cathode = 0.032
        self.w = 0.224
        self.p = 0.512
        self.ysize = 0

        pos_1_box1 = self.pos_0 - np.array([0, 0, self.h])
        pos_2_box1 = self.pos_0 + np.array([self.w, self.ysize, -self.h_cathode-self.h])

        pos_1_box2 = self.pos_0 + np.array([self.w, 0, 0])
        pos_2_box2 = self.pos_0 + np.array([self.p, self.ysize, -self.h_cathode])
        self.geom = Geometry.Rectangular(pos_1_box1, pos_2_box1, substrate_normales_indx = [2])

    def test_is_outside(self):

        self.assertEqual(self.geom._is_outside(np.array([0, 0, 0])), Geometry.STATUS['Exit'])
        self.assertEqual(self.geom._is_outside(self.pos_0 + np.array([self.w/2, 0, -self.h-0.005])), Geometry.STATUS['Inside'])
        self.assertEqual(self.geom._is_outside(self.pos_0 + np.array([self.w/2, 500, -self.h-0.005])), Geometry.STATUS['Inside'])
        self.assertEqual(self.geom._is_outside(self.pos_0 + np.array([2*self.w, 500, -self.h-0.005])), Geometry.STATUS['Exit'])
        self.assertEqual(self.geom._is_outside(np.array([0, 500, 0])), Geometry.STATUS['Exit'])

    def test_get_outer_way(self):

        self.assertEqual(self.geom.get_outer_way(np.array([0, 0, 0])), 2)
        self.assertEqual(self.geom.get_outer_way(np.array([0, 4, 0])), 2)

        self.assertEqual(self.geom.get_outer_way(np.array([0.010, 0, 0.512 - self.h - 2*self.h_cathode])), 2)
        self.assertEqual(self.geom.get_outer_way(np.array([0.010, 600, 0.512 - self.h - 2*self.h_cathode])), 2)

        self.assertEqual(self.geom.get_outer_way(np.array([self.w + 0.005, 0, 0.512 - self.h - self.h_cathode + 0.005])), 1)
        self.assertEqual(self.geom.get_outer_way(np.array([self.w + 0.005, 100, 0.512 - self.h - self.h_cathode + 0.005])), 1)

        self.assertEqual(self.geom.get_outer_way(np.array([-0.005, 0, 0.512 - self.h - self.h_cathode + 0.005])), 0)
        self.assertEqual(self.geom.get_outer_way(np.array([-0.005, -2, 0.512 - self.h - self.h_cathode + 0.005])), 0)

        self.assertEqual(self.geom.get_outer_way(np.array([0.01, 0, 0.512 - self.h + 0.005])), 3)
        self.assertEqual(self.geom.get_outer_way(np.array([0.01, 200, 0.512 - self.h + 0.005])), 3)
        

if __name__ == '__main__':
    unittest.main()
