import unittest
import electron
import ElTransport
import numpy as np

def _initial_process_single_electron(coor, energy, effective_mass):

    single_electron = electron.Electrons(coor[0], coor[1], coor[2], 0, 0, 0)
    single_electron.set_electron_properties(effective_mass)

    veloicity = ElTransport.make_initial_dir(single_electron, energy)
    single_electron.set_veloicity(veloicity)

    return single_electron

class Test_electron(unittest.TestCase):

    def setUp(self):

        coor = np.zeros(3)
        self.electron1 = _initial_process_single_electron(coor, 1, 0.1)
        self.electron2 = _initial_process_single_electron(coor, 0, 0.1)
        self.electron3 = _initial_process_single_electron(coor, 2, 0.1)

    def test_get_E(self):

        self.assertEqual(round(self.electron1.get_E(), 2), 1)
        self.assertEqual(round(self.electron2.get_E(), 2), 0)
        self.assertEqual(round(self.electron3.get_E(), 2), 2)
    
    def test_add_energy(self):

        e01 = 1
        e02 = 0.1
        e03 = 2

        coor = np.zeros(3)
        electron1 = _initial_process_single_electron(coor, e01, 0.1)
        electron2 = _initial_process_single_electron(coor, e02, 0.1)
        electron3 = _initial_process_single_electron(coor, e03, 0.1)

        deltaE1 = -0.1
        deltaE2 = 0.1
        deltaE3 = 0.1

        electron1.add_energy(deltaE1)
        electron2.add_energy(deltaE2)
        electron3.add_energy(deltaE3)

        self.assertEqual(round(electron1.get_E(), 2), e01 + deltaE1)
        self.assertEqual(round(electron2.get_E(), 2), e02 + deltaE2)
        self.assertEqual(round(electron3.get_E(), 2), e03 + deltaE3)

if __name__ == '__main__':
    unittest.main()
