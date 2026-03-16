import numpy as np
import pytest
from unittest.mock import Mock, patch
from electron import Electrons
import ElTransport

class TestElTransport:
    
    def test_make_new_coor(self):
        """Test updating electron coordinates"""
        electrons = Electrons(3)
        electrons.set_velocity(np.array([[1.0, 2.0, 3.0],
                                         [4.0, 5.0, 6.0],
                                         [7.0, 8.0, 9.0]]))
        
        electrons.set_xcoor(np.array([[0.0, 0.0, 0.0],
                                      [1.0, 1.0, 1.0],
                                      [2.0, 2.0, 2.0]]))
        
        electrons.set_flags(np.array([1, 1, 1]))
        
        dt = 0.1
        ElTransport._make_new_coor(electrons, dt)
        
        expected = np.array([[0.1, 0.2, 0.3],
                             [1.4, 1.5, 1.6],
                             [2.7, 2.8, 2.9]])
        np.testing.assert_array_almost_equal(electrons.get_xcoor(), expected)
    
    def test_make_new_coor_with_inactive_electrons(self):
        """Test: inactive electrons don't move"""
        electrons = Electrons(2)
        electrons.set_velocity(np.array([[1.0, 1.0, 1.0],
                                         [2.0, 2.0, 2.0]]))
        electrons.set_flags(np.array([0, 1]))  # First is inactive
        
        dt = 1.0
        ElTransport._make_new_coor(electrons, dt)
        
        coords = electrons.get_xcoor()
        assert coords[0, 0] == 0.0  # First electron
        assert coords[1, 0] == 2.0  # Second electron
    
    def test_make_p_l_e_e(self):
        """Test scattering probability calculation"""
        electrons = Electrons(2)
        electrons.set_electron_properties(0.1)
        electrons.set_velocity(np.array([[1.0, 0.0, 0.0],
                                         [2.0, 0.0, 0.0]]))
        
        def mock_l_e_e(energy):
            return 1.0  # constant mean free path
        
        dt = 0.5
        p = ElTransport._make_p_l_e_e(electrons, mock_l_e_e, dt)
        
        v = electrons.get_module_velocity()
        expected = 1 - np.exp(-v * dt / 1.0)
        np.testing.assert_array_almost_equal(p, expected)
    
    def test_is_scat(self):
        """Test scattering occurrence check"""
        with patch('numpy.random.rand', return_value=0.3):
            assert ElTransport._is_scat(0.5) == True
            assert ElTransport._is_scat(0.3) == False
            assert ElTransport._is_scat(0.1) == False
    
    def test_make_velocities_after_scattering(self):
        """Test velocity update after scattering"""
        electrons = Electrons(2)
        electrons.set_electron_properties(0.1)
        
        initial_velocity = np.array([[1.0, 0.0, 0.0],
                                     [0.0, 2.0, 0.0]])
        electrons.set_velocity(initial_velocity)
        
        initial_energy = electrons.get_E()
        
        delta_E = np.array([[1.0], [2.0]])
        
        any_scattering = np.array([True, True])
        N_el_have_scat = 2
        
        np.random.seed(42)
        ElTransport._make_velocities_after_scattering(
            electrons, delta_E, any_scattering, N_el_have_scat
        )
        
        new_energy = electrons.get_E()
        assert np.all(new_energy > initial_energy)
        
        module_vel = electrons.get_module_velocity()
        expected_module_vel = 1e-3 * np.sqrt(
            (initial_energy.reshape(-1, 1) + delta_E) * 2 * ElTransport.EV_CONST / 
            (0.1 * ElTransport.M_E)
        ).flatten()
        
        np.testing.assert_array_almost_equal(module_vel, expected_module_vel, decimal=5)
    
    def test_make_scatterings(self):
        """Test scattering process"""
        electrons = Electrons(3)
        electrons.set_electron_properties(0.1)
        electrons.set_flags(np.array([1, 1, 1]))
        
        electrons.set_velocity(np.array([[1.0, 0.0, 0.0],
                                         [2.0, 0.0, 0.0],
                                         [3.0, 0.0, 0.0]]))
        
        scatterings_l_e_e = [
            lambda energy: 1.0,
            lambda energy: 2.0
        ]
        
        scatterings_E = [0.1, 0.2]
        
        dt = 0.5
        
        with patch('ElTransport._make_p_l_e_e') as mock_p, \
             patch('numpy.random.rand') as mock_rand:
            
            mock_p.return_value = np.array([0.8, 0.8, 0.8])
            mock_rand.return_value = 0.1
            
            with patch('ElTransport._make_velocities_after_scattering') as mock_make_vel:
                ElTransport._make_scatterings(electrons, dt, scatterings_l_e_e, scatterings_E)
                
                mock_make_vel.assert_called_once()
                
                args, kwargs = mock_make_vel.call_args
                assert args[0] is electrons
                assert args[1].shape == (3, 1)
                assert np.all(args[2] == [True, True, True])
                assert args[3] == 3
    
    def test_make_scatterings_no_scattering(self):
        """Test: when no scattering occurs"""
        electrons = Electrons(2)
        electrons.set_flags(np.array([1, 1]))
        
        scatterings_l_e_e = [lambda energy: 1.0]
        scatterings_E = [0.1]
        
        dt = 0.5
        
        with patch('ElTransport._make_p_l_e_e') as mock_p, \
             patch('numpy.random.rand') as mock_rand:
            
            mock_p.return_value = np.array([0.1, 0.1])
            mock_rand.return_value = 0.5
            
            with patch('ElTransport._make_velocities_after_scattering') as mock_make_vel:
                ElTransport._make_scatterings(electrons, dt, scatterings_l_e_e, scatterings_E)
                
                mock_make_vel.assert_not_called()
    
    def test_make_initial_velocity(self):
        """Test initial velocity creation"""
        electrons = Electrons(100)
        electrons.set_electron_properties(0.2)
        
        energy = 1.0
        
        np.random.seed(42)
        ElTransport.make_initial_veloicity(electrons, energy)
        
        velocities = electrons.get_velocity()
        assert velocities.shape == (100, 3)
        
        module_vel = electrons.get_module_velocity()
        expected_module_vel = 1e-3 * np.sqrt(
            energy * 2 * ElTransport.EV_CONST / 
            (0.2 * ElTransport.M_E)
        )
        
        np.testing.assert_allclose(module_vel.mean(), expected_module_vel, rtol=0.1)
        
        norms = np.linalg.norm(velocities / module_vel.reshape(-1, 1), axis=1)
        np.testing.assert_array_almost_equal(norms, 1.0)
    
    def test_transport_process(self):
        """Integration test of transport process"""
        electrons = Electrons(10)
        electrons.set_electron_properties(0.15)
        electrons.set_flags(np.ones(10))
        
        scatterings_l_e_e = [lambda energy: 1.0]
        scatterings_E = [0.1]
        
        dt = 0.01
        
        with patch('ElTransport._make_new_coor') as mock_new_coor, \
             patch('ElTransport._make_scatterings') as mock_scatterings:
            
            ElTransport.transport_process(electrons, dt, scatterings_l_e_e, scatterings_E)
            
            mock_new_coor.assert_called_once_with(electrons, dt)
            mock_scatterings.assert_called_once_with(electrons, dt, scatterings_l_e_e, scatterings_E)
    
    def test_constants(self):
        """Test constants"""
        assert ElTransport.C_CONST == 2.99792458
        assert ElTransport.EV_CONST == 1.602176634
        assert ElTransport.M_E == 9.109
    
    def test_edge_cases(self):
        """Test edge cases"""
        # Zero time
        electrons = Electrons(2)
        electrons.set_velocity(np.array([[1.0, 0.0, 0.0],
                                         [0.0, 1.0, 0.0]]))
        
        initial_coords = electrons.get_xcoor().copy()
        dt = 0.0
        
        ElTransport._make_new_coor(electrons, dt)
        
        np.testing.assert_array_equal(electrons.get_xcoor(), initial_coords)
        
        # Zero velocity
        electrons.set_velocity(np.array([[0.0, 0.0, 0.0],
                                         [0.0, 0.0, 0.0]]))
        electrons.set_xcoor(initial_coords)
        
        ElTransport._make_new_coor(electrons, 1.0)
        
        np.testing.assert_array_equal(electrons.get_xcoor(), initial_coords)
    
    def test_velocity_calculation_with_energy(self):
        """Test velocity calculation from energy"""
        electrons = Electrons(1)
        electrons.set_electron_properties(0.1)
        
        # Set velocity manually
        test_velocity = np.array([[1000.0, 0.0, 0.0]])  # 1e3 m/s
        electrons.set_velocity(test_velocity)
        
        # Calculate energy
        energy = electrons.get_E()
        
        # Now set new velocity using make_initial_veloicity with this energy
        np.random.seed(42)
        ElTransport.make_initial_veloicity(electrons, energy[0])
        
        # Check if energy is approximately conserved
        new_energy = electrons.get_E()
        np.testing.assert_allclose(new_energy[0], energy[0], rtol=0.1)

if __name__ == "__main__":
    pytest.main([__file__, "-v"])