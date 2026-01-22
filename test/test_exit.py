import unittest
import Geometry
import numpy as np
import electron
import ElectronExit as elex
import MyScatterings as scat
import ElTransport as eltrans
import sys
import os
import pytest
import numpy as np
from unittest.mock import Mock, patch
import ElectronExit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import the module to test
import ElectronExit

class TestPEexit:
    """Test the _p_exit function"""
    
    def test_p_exit_valid_conditions(self):
        """Test probability calculation when exit conditions are met"""
        E = 10.0
        E_a = 2.0
        cos_angle = 0.8
        
        result = ElectronExit._p_exit(E, E_a, cos_angle)
        
        E_exit = E * cos_angle * cos_angle
        expected = 4 * np.sqrt(E_exit * np.abs(E_exit - E_a)) / (np.sqrt(np.abs(E_exit - E_a)) + np.sqrt(E_exit))**2
        
        assert np.isclose(result, expected)
        assert result > 0
    
    def test_p_exit_invalid_conditions(self):
        """Test probability calculation when exit conditions are not met"""
        E = 10.0
        E_a = 8.0
        cos_angle = 0.2  # Too small cosine angle
        
        result = ElectronExit._p_exit(E, E_a, cos_angle)
        
        assert result == 0
    
    def test_p_exit_edge_case(self):
        """Test probability calculation at edge conditions"""
        E = 10.0
        E_a = 5.0
        cos_angle = np.sqrt(E_a / E)  # Exactly at threshold
        
        result = ElectronExit._p_exit(E, E_a, cos_angle)
        
        # Should be 0 because condition is sqrt(E_a/E) < cos_angle (not <=)
        assert result == 0
    
    def test_p_exit_array_input(self):
        """Test probability calculation with array inputs"""
        E = np.array([10.0, 20.0, 30.0])
        E_a = 5.0
        cos_angle = np.array([0.9, 0.6, 0.3])
        
        result = ElectronExit._p_exit(E, E_a, cos_angle)
        
        assert isinstance(result, np.ndarray)
        assert len(result) == 3
        assert result[0] > 0  # First should exit
        assert result[2] == 0  # Third should not exit

class TestIsExit:
    """Test the is_exit function"""
    
    def test_is_exit_true(self):
        """Test when electron should exit"""
        mock_geom = Mock()
        mock_electrons = Mock()
        mock_semiconductor = Mock()
        
        mock_electrons.get_E.return_value = 10.0
        mock_semiconductor.get_E_a.return_value = 2.0
        mock_geom.get_cos_angle.return_value = 0.9
        
        # Probability will be ~0.847, so with random 0.5 it should exit
        with patch('numpy.random.rand', return_value=0.5):
            result = ElectronExit.is_exit(mock_geom, mock_electrons, mock_semiconductor, 0, 0, True)
            
            assert result == True
            mock_electrons.get_E.assert_called_once_with(0)
            mock_semiconductor.get_E_a.assert_called_once()
            mock_geom.get_cos_angle.assert_called_once_with(mock_electrons, 0, 0)
    
    def test_is_exit_false(self):
        """Test when electron should not exit"""
        mock_geom = Mock()
        mock_electrons = Mock()
        mock_semiconductor = Mock()
        
        mock_electrons.get_E.return_value = 2.1
        mock_semiconductor.get_E_a.return_value = 2.0
        mock_geom.get_cos_angle.return_value = 0.85
        
        # Probability will be ~0.847, so with random 0.9 it should not exit
        with patch('numpy.random.rand', return_value=0.9):
            result = ElectronExit.is_exit(mock_geom, mock_electrons, mock_semiconductor, 0, 0, True)
            
            assert result == False

class TestGetThisElectronEmittance:
    """Test the _get_this_electron_emitance function"""
    
    def test_get_this_electron_emittance(self):
        """Test emittance calculation for a single electron"""
        mock_geom = Mock()
        mock_electrons = Mock()
        
        mock_electrons.get_E.return_value = 10.0
        mock_geom.get_cos_angle.return_value = 0.6
        
        result = ElectronExit._get_this_electron_emitance(mock_geom, mock_electrons, 0, 0)
        
        expected = 10.0 * (1 - 0.6 * 0.6)
        assert np.isclose(result, expected)
        
        mock_electrons.get_E.assert_called_once_with(0)
        mock_geom.get_cos_angle.assert_called_once_with(mock_electrons, 0, 0)

class TestReflectionProcess:
    """Test the _reflection_process function"""
    
    def test_reflection_process(self):
        """Test reflection process updates electron coordinates"""
        mock_geom = Mock()
        mock_electrons = Mock()
        
        new_coor = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        mock_geom.get_new_coors_after_reflect.return_value = new_coor
        
        ElectronExit._reflection_process(mock_geom, mock_electrons, 0, 0)
        
        mock_geom.get_new_coors_after_reflect.assert_called_once_with(mock_electrons, 0, 0)
        mock_electrons.set_coor.assert_called_once_with(new_coor, 0)

class TestExitProcess:
    """Test the exit_process function"""
    
    def test_exit_process_no_electrons(self):
        """Test exit process with no electrons alive"""
        mock_geom = Mock()
        mock_electrons = Mock()
        mock_semiconductor = Mock()
        
        mock_electrons.get_N_el.return_value = 3
        mock_electrons.is_alive.side_effect = [False, False, False]
        
        result = ElectronExit.exit_process(mock_geom, mock_electrons, mock_semiconductor)
        
        assert result['N_exit'] == 0
        assert result['Emmitance'] == 0
        assert mock_electrons.get_N_el.call_count == 1
        assert mock_electrons.is_alive.call_count == 3
    
    def test_exit_process_electron_exits(self):
        """Test exit process when an electron exits"""
        mock_geom = Mock()
        mock_electrons = Mock()
        mock_semiconductor = Mock()
        
        mock_electrons.get_N_el.return_value = 2
        mock_electrons.is_alive.side_effect = [True, False]
        mock_geom.get_out_plane.return_value = 0
        mock_geom.get_type_exit_plane.return_value = 1
        
        # Mock the is_exit function to return True
        with patch.object(ElectronExit, 'is_exit', return_value=True):
            with patch.object(ElectronExit, '_get_this_electron_emitance', return_value=5.0):
                result = ElectronExit.exit_process(mock_geom, mock_electrons, mock_semiconductor)
        
        assert result['N_exit'] == 1
        assert result['Emmitance'] == 5.0
        mock_electrons.kill_electron.assert_called_once_with(0)
    
    def test_exit_process_electron_reflects(self):
        """Test exit process when an electron reflects"""
        mock_geom = Mock()
        mock_electrons = Mock()
        mock_semiconductor = Mock()
        
        mock_electrons.get_N_el.return_value = 2
        mock_electrons.is_alive.side_effect = [True, False]
        mock_geom.get_out_plane.return_value = 0
        mock_geom.get_type_exit_plane.return_value = 1
        
        # Mock the is_exit function to return False (reflection)
        with patch.object(ElectronExit, 'is_exit', return_value=False):
            with patch.object(ElectronExit, '_reflection_process') as mock_reflection:
                result = ElectronExit.exit_process(mock_geom, mock_electrons, mock_semiconductor)
        
        assert result['N_exit'] == 0
        assert result['Emmitance'] == 0
        mock_reflection.assert_called_once_with(mock_geom, mock_electrons, 0, 0)
    
    def test_exit_process_electron_dies(self):
        """Test exit process when an electron dies on non-exit plane"""
        mock_geom = Mock()
        mock_electrons = Mock()
        mock_semiconductor = Mock()
        
        mock_electrons.get_N_el.return_value = 2
        mock_electrons.is_alive.side_effect = [True, False]
        mock_geom.get_out_plane.return_value = 0
        mock_geom.get_type_exit_plane.return_value = 0  # Non-exit plane
        
        result = ElectronExit.exit_process(mock_geom, mock_electrons, mock_semiconductor)
        
        assert result['N_exit'] == 0
        assert result['Emmitance'] == 0
        mock_electrons.kill_electron.assert_called_once_with(0)
    
    def test_exit_process_electron_inside(self):
        """Test exit process when electron is inside (no out plane)"""
        mock_geom = Mock()
        mock_electrons = Mock()
        mock_semiconductor = Mock()
        
        mock_electrons.get_N_el.return_value = 2
        mock_electrons.is_alive.side_effect = [True, False]
        mock_geom.get_out_plane.return_value = -1  # Inside, no exit
        
        result = ElectronExit.exit_process(mock_geom, mock_electrons, mock_semiconductor)
        
        assert result['N_exit'] == 0
        assert result['Emmitance'] == 0
        mock_geom.get_type_exit_plane.assert_not_called()
        mock_electrons.kill_electron.assert_not_called()
    
    def test_exit_process_multiple_electrons(self):
        """Test exit process with multiple electrons in different states"""
        mock_geom = Mock()
        mock_electrons = Mock()
        mock_semiconductor = Mock()
        
        mock_electrons.get_N_el.return_value = 4
        mock_electrons.is_alive.side_effect = [True, False, True, True]
        mock_geom.get_out_plane.side_effect = [0, -1, 1, 0]  # Different out planes
        mock_geom.get_type_exit_plane.return_value = 1
        
        # Mock is_exit to return True for first, False for third and fourth
        with patch.object(ElectronExit, 'is_exit', side_effect=[True, False, False]):
            with patch.object(ElectronExit, '_get_this_electron_emitance', return_value=3.0):
                with patch.object(ElectronExit, '_reflection_process'):
                    result = ElectronExit.exit_process(mock_geom, mock_electrons, mock_semiconductor)
        
        assert result['N_exit'] == 1
        assert result['Emmitance'] == 3.0
        assert mock_electrons.kill_electron.call_count == 1

if __name__ == '__main__':
    pytest.main([__file__, '-v'])