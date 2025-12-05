import numpy as np
import pytest
import unittest
from electron import Electrons, _check_sizes


class TestCheckSizesFunction(unittest.TestCase):
    def test_check_sizes_valid(self):
        N_el = 3
        mass2 = np.zeros((3, 2))
        _check_sizes(N_el, mass2)
    
    def test_check_sizes_invalid(self):
        N_el = 3
        mass2 = np.zeros((4, 2))
        with pytest.raises(ValueError, match='Not equal sizes'):
            _check_sizes(N_el, mass2)


class TestElectronsInitialization(unittest.TestCase):
    def test_init_valid(self):
        N = 5
        electrons = Electrons(N)
        
        self.assertEqual(electrons.N_el_in_array, N)
        self.assertEqual(electrons.coor.shape, (N, 7))
        np.testing.assert_array_equal(electrons.coor, np.zeros((N, 7)))
    
    def test_init_zero_or_negative(self):
        with pytest.raises(ValueError, match='N_el_in_array must be greater than zero'):
            Electrons(0)
        
        with pytest.raises(ValueError, match='N_el_in_array must be greater than zero'):
            Electrons(-5)


class TestElectronsProperties(unittest.TestCase):
    
    def setUp(self):
        self.N = 3
        self.electrons = Electrons(self.N)
    
    def test_set_electron_properties_valid(self):
        mass = 0.5
        self.electrons.set_electron_properties(mass)
        self.assertEqual(self.electrons.effective_mass, mass)
    
    def test_set_electron_properties_invalid(self):
        with pytest.raises(ValueError, match='Effective mass must be positive'):
            self.electrons.set_electron_properties(0)
        
        with pytest.raises(ValueError, match='Effective mass must be positive'):
            self.electrons.set_electron_properties(-1.0)


class TestElectronsCoordinates(unittest.TestCase):
    
    def setUp(self):
        self.N = 3
        self.electrons = Electrons(self.N)
        test_coor = np.array([
            [1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0],
            [4.0, 5.0, 6.0, 0.4, 0.5, 0.6, 0],
            [7.0, 8.0, 9.0, 0.7, 0.8, 0.9, 0]
        ])
        self.electrons.coor = test_coor.copy()
    
    def test_set_coor_all(self):
        new_coor = np.array([
            [10.0, 11.0, 12.0, 0.1, 0.2, 0.3],
            [13.0, 14.0, 15.0, 0.4, 0.5, 0.6],
            [16.0, 17.0, 18.0, 0.7, 0.8, 0.9]
        ])
        self.electrons.set_coor(new_coor)
        
        expected = np.hstack([new_coor, np.zeros((self.N, 1))])
        np.testing.assert_array_almost_equal(self.electrons.coor, expected)
    
    def test_set_coor_single(self):
        new_coor = np.array([10.0, 11.0, 12.0, 0.1, 0.2, 0.3])
        self.electrons.set_coor(new_coor, indx_el=1)
        
        np.testing.assert_array_almost_equal(self.electrons.coor[1, :-1], new_coor)
        np.testing.assert_array_almost_equal(self.electrons.coor[0, :-1], 
                                           [1.0, 2.0, 3.0, 0.1, 0.2, 0.3])
    
    def test_set_coor_wrong_size(self):
        wrong_coor = np.array([
            [10.0, 11.0, 12.0, 0.1, 0.2, 0.3],
            [13.0, 14.0, 15.0, 0.4, 0.5, 0.6]
        ])  
        
        with pytest.raises(ValueError, match='Not equal sizes'):
            self.electrons.set_coor(wrong_coor)
    
    def test_get_coor_all(self):
        result = self.electrons.get_coor()
        expected = self.electrons.coor[:, :-1]
        np.testing.assert_array_almost_equal(result, expected)
    
    def test_get_coor_single(self):
        result = self.electrons.get_coor(indx_el=1)
        expected = self.electrons.coor[1, :-1]
        np.testing.assert_array_almost_equal(result, expected)
    
    def test_set_xcoor_all(self):
        new_xcoor = np.array([
            [10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0],
            [16.0, 17.0, 18.0]
        ])
        self.electrons.set_xcoor(new_xcoor)
        np.testing.assert_array_almost_equal(self.electrons.coor[:, :3], new_xcoor)
        np.testing.assert_array_almost_equal(self.electrons.coor[:, 3:-1], 
                                           self.electrons.coor[:, 3:-1])
    
    def test_add_coor(self):
        shift = np.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0]
        ])
        original_coor = self.electrons.coor[:, :3].copy()
        self.electrons.add_coor(shift)
        
        expected = original_coor + shift
        np.testing.assert_array_almost_equal(self.electrons.coor[:, :3], expected)


class TestElectronsVelocity(unittest.TestCase):
    
    def setUp(self):
        self.N = 3
        self.electrons = Electrons(self.N)
        test_coor = np.array([
            [1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0],
            [4.0, 5.0, 6.0, 0.4, 0.5, 0.6, 0],
            [7.0, 8.0, 9.0, 0.7, 0.8, 0.9, 0]
        ])
        self.electrons.coor = test_coor.copy()
        self.electrons.set_electron_properties(0.5)
    
    def test_set_velocity_all(self):
        new_velocity = np.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0]
        ])
        self.electrons.set_velocity(new_velocity)
        
        expected = np.hstack([
            self.electrons.coor[:, :3],
            new_velocity,
            self.electrons.coor[:, -1:]
        ])
        np.testing.assert_array_almost_equal(self.electrons.coor, expected)
    
    def test_get_velocity(self):
        result_all = self.electrons.get_velocity()
        expected_all = self.electrons.coor[:, 3:-1]
        np.testing.assert_array_almost_equal(result_all, expected_all)
        
        result_single = self.electrons.get_velocity(indx_el=1)
        expected_single = self.electrons.coor[1, 3:-1]
        np.testing.assert_array_almost_equal(result_single, expected_single)
    
    def test_get_module_velocity(self):

        module_all = self.electrons.get_module_velocity()
        expected_all = np.sqrt(0.1**2 + 0.2**2 + 0.3**2)
        self.assertAlmostEqual(module_all[0], expected_all, places=10)

        module_single = self.electrons.get_module_velocity(indx_el=2)
        expected_single = np.sqrt(0.7**2 + 0.8**2 + 0.9**2)
        self.assertAlmostEqual(module_single, expected_single, places=10)


class TestElectronsEnergy(unittest.TestCase):
    
    def setUp(self):
        self.N = 2
        self.electrons = Electrons(self.N)

        test_coor = np.array([
            [0.0, 0.0, 0.0, 1e6, 0.0, 0.0, 0], 
            [0.0, 0.0, 0.0, 0.0, 2e6, 0.0, 0]  
        ])
        self.electrons.coor = test_coor.copy()
        self.electrons.set_electron_properties(0.5)
    
    def test_get_energy(self):

        energy1 = self.electrons.get_E(indx_el=0)

        expected1 = 0.5 * 9.109 * (1e6)**2 / (2 * 1.602176634) * 1e6
        self.assertAlmostEqual(energy1, expected1, delta=expected1*0.01)
        
        energies = self.electrons.get_E()
        self.assertEqual(len(energies), self.N)
    
    def test_add_energy_all(self):

        initial_energy = self.electrons.get_E(0)
        delta_E = 1e-18
        
        self.electrons.add_energy(np.array([delta_E, delta_E]))
        
        final_energy = self.electrons.get_E(0)


        self.assertAlmostEqual(final_energy, initial_energy + delta_E, 
                             delta=final_energy*0.01)
    
    def test_add_energy_single(self):

        initial_energy = self.electrons.get_E(indx_el=1).copy()
        delta_E = 1e-18
        
        self.electrons.add_energy(delta_E, indx_el=1)
        
        final_energy = self.electrons.get_E(indx_el=1)
        self.assertAlmostEqual(final_energy, initial_energy + delta_E,
                             delta=final_energy*0.01)


class TestElectronsFlags(unittest.TestCase):
    
    def setUp(self):
        self.N = 4
        self.electrons = Electrons(self.N)
    
    def test_set_get_flags(self):

        flags = np.array([0, 1, 0, 1])
        self.electrons.set_flags(flags)
        np.testing.assert_array_equal(self.electrons.get_flags(), flags)

        self.electrons.set_flags(1, indx_el=2)
        self.assertEqual(self.electrons.get_flags(indx_el=2), 1)
    
    def test_kill_electron(self):
        self.electrons.kill_electron(2)
        self.assertEqual(self.electrons.coor[2, -1], 1)
        self.assertFalse(self.electrons.is_alive(2))
    
    def test_is_alive(self):

        self.electrons.set_flags(0, indx_el=0)
        self.electrons.set_flags(1, indx_el=1)
        
        self.assertTrue(self.electrons.is_alive(0))
        self.assertFalse(self.electrons.is_alive(1))
    
    def test_is_end(self):

        self.electrons.set_flags(np.array([0, 1, 0, 1]))
        self.assertFalse(self.electrons.is_end())
        
        self.electrons.set_flags(np.array([1, 1, 1, 1]))
        self.assertTrue(self.electrons.is_end())
    
    def test_kill_low_energy_electron(self):

        self.N = 3
        self.electrons = Electrons(self.N)
        self.electrons.set_electron_properties(0.5)
        
        velocities = np.array([
            [1e4, 0.0, 0.0],  
            [1e6, 0.0, 0.0], 
            [1e4, 0.0, 0.0]  
        ])
        self.electrons.set_velocity(velocities)
        
        low_energy = self.electrons.get_E(indx_el=0)
        high_energy = self.electrons.get_E(indx_el=1)
        kill_energy = (low_energy + high_energy) / 2
        
        self.electrons.kill_low_energy_electron(kill_energy)
        
        flags = self.electrons.get_flags()
        expected = np.array([1, 0, 1]) 
        np.testing.assert_array_equal(flags, expected)


class TestElectronsOtherMethods(unittest.TestCase):
    
    def setUp(self):
        self.N = 2
        self.electrons = Electrons(self.N)
        self.electrons.set_electron_properties(0.5)
    
    def test_get_N_el_in_ar(self):
        self.assertEqual(self.electrons.get_N_el_in_ar(), self.N)
    
    def test_get_effective_mass(self):
        self.assertEqual(self.electrons.get_effective_mass(), 0.5)
    
    def test_get_log(self):
        self.electrons.E_a = 0.1
        self.electrons.E_g = 1.5
        self.electrons.delta_E_DOS = 0.05
        
        log = self.electrons.get_log()
        expected = f"M_e={0.5}\nE_a={0.1}\nE_g={1.5}\ndelta_E_DOS={0.05}\n"
        self.assertEqual(log, expected)


if __name__ == '__main__':
    unittest.main()