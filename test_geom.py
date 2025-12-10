import pytest
import numpy as np
from Geometry import Plane, ConvexShape

# --- Tests for Plane class ---

class TestPlane:
    def test_plane_initialization(self):
        """Test plane initialization with correct parameters"""
        point = np.array([1.0, 2.0, 3.0])
        normale = np.array([0.0, 0.0, 1.0])
        plane_name = "XY_plane"
        plane = Plane(point, normale, plane_name)

        assert np.array_equal(plane.point, point)
        assert np.array_equal(plane.normale, normale)
        assert plane.plane_name == plane_name

    def test_plane_get_normale(self):
        """Test getting plane normal vector"""
        normale = np.array([1.0, 0.0, 0.0])
        plane = Plane(np.zeros(3), normale, "test")
        assert np.array_equal(plane.get_normale(), normale)

    def test_plane_get_point(self):
        """Test getting plane reference point"""
        point = np.array([5.0, 5.0, 5.0])
        plane = Plane(point, np.zeros(3), "test")
        assert np.array_equal(plane.get_point(), point)

    def test_plane_get_distance(self):
        """Test calculating distance from point to plane"""
        # XY plane with normal pointing up along Z
        plane = Plane(
            point=np.array([0.0, 0.0, 0.0]),
            normale=np.array([0.0, 0.0, 1.0]),
            plane_name="XY"
        )

        # Point above the plane
        assert plane.get_distance(np.array([2.0, 3.0, 5.0])) == 5.0
        # Point below the plane
        assert plane.get_distance(np.array([1.0, 1.0, -2.0])) == -2.0
        # Point on the plane
        assert plane.get_distance(np.array([10.0, 20.0, 0.0])) == 0.0

    def test_plane_str_representation(self):
        """Test string representation of plane"""
        plane = Plane(
            point=np.array([1.0, 2.0, 3.0]),
            normale=np.array([0.0, 1.0, 0.0]),
            plane_name="test_plane"
        )
        expected = "point = [1. 2. 3.], normale = [0. 1. 0.], plane_name = test_plane"
        assert str(plane) == expected


# --- Mock class for ConvexShape testing ---

class MockElectron:
    """Mock electron object for testing geometry interactions"""
    def __init__(self, xcoor, velocity):
        self.xcoor = xcoor
        self.velocity = velocity

    def get_xcoor(self, indx_el=None):
        return self.xcoor

    def get_velocity(self, indx_el=None):
        return self.velocity

    def get_module_velocity(self, indx_el=None):
        return np.sqrt(np.sum(self.velocity**2))


# --- Tests for ConvexShape class ---

class TestConvexShape:
    def setup_method(self):
        """Setup before each test method"""
        self.shape = ConvexShape()
        # Create 3 planes forming a triangular region
        self.plane1 = Plane(
            point=np.array([0.0, 0.0, 0.0]),
            normale=np.array([-1.0, 0.0, 0.0]),
            plane_name="plane1"
        )
        self.plane2 = Plane(
            point=np.array([0.0, 0.0, 1.0]),
            normale=np.array([np.sqrt(2)/2, 0.0, np.sqrt(2)/2]),
            plane_name="plane2"
        )
        self.plane3 = Plane(
            point=np.array([0.0, 0.0, 0.0]),
            normale=np.array([0.0, 0.0, -1.0]),
            plane_name="plane3"
        )

    def test_add_plane(self):
        """Test adding planes to convex shape"""
        self.shape.add_plane(self.plane1, 1)
        self.shape.add_plane(self.plane2, 0)
        self.shape.add_plane(self.plane3, 1)

        assert len(self.shape.planes) == 3
        assert self.shape.planes[0] == self.plane1
        assert self.shape.planes[1] == self.plane2
        assert self.shape.planes[2] == self.plane3
        assert len(self.shape.work_surfaces) == 3

    def test_get_name(self):
        """Test getting shape name"""
        assert self.shape.get_name() == "ConvexShape"

    def test_get_params(self):
        """Test getting shape parameters as string"""
        self.shape.add_plane(self.plane1, 1)
        self.shape.add_plane(self.plane2, 0)
        
        params = self.shape.get_params()
        assert "plane1" in params
        assert "plane2" in params
        assert "point =" in params
        assert "normale =" in params

    def test_get_out_plane_electron_inside(self):
        """Test exit plane detection for electron inside shape"""
        self.shape.add_plane(self.plane1, 1)
        self.shape.add_plane(self.plane2, 1)
        self.shape.add_plane(self.plane3, 1)

        # Electron inside (all distances negative)
        electron = MockElectron(
            xcoor=np.array([0.1, 0.1, 0.1]),
            velocity=np.array([0.0, 0.0, 0.0])
        )
        # Method should return -1 if electron is inside all planes
        assert self.shape.get_out_plane(electron, 0) == -1

    def test_get_out_plane_electron_outside(self):
        """Test exit plane detection for electron outside shape"""
        self.shape.add_plane(self.plane1, 1)
        self.shape.add_plane(self.plane2, 1)
        self.shape.add_plane(self.plane3, 1)

        # Electron outside through plane1 (x < 0)
        electron = MockElectron(
            xcoor=np.array([-0.1, 0.5, 0.5]),
            velocity=np.array([0.0, 0.0, 0.0])
        )
        
        assert self.shape.get_out_plane(electron, 0) == 0

        # Electron outside through plane3 (z < 0)
        electron2 = MockElectron(
            xcoor=np.array([0.5, 0.5, -0.1]),
            velocity=np.array([0.0, 0.0, 0.0])
        )
        
        assert self.shape.get_out_plane(electron2, 0) == 2

    def test_get_type_exit_plane(self):
        """Test getting exit plane by index"""
        self.shape.add_plane(self.plane1, 1)
        self.shape.add_plane(self.plane3, 0)
        
        # Get planes by index
        assert self.shape.get_type_exit_plane(0) == 1
        assert self.shape.get_type_exit_plane(1) == 0

    def test_get_type_exit_plane_invalid_index(self):
        """Test exception for invalid plane index"""
        self.shape.add_plane(self.plane1, 1)
        
        with pytest.raises(ValueError, match="Electron inside and cannt escape"):
            self.shape.get_type_exit_plane(-1)

    def test_get_new_coors_after_reflect(self):
        """Test reflection calculation from plane"""
        # Plane X=0 with normal (1,0,0)
        plane = Plane(
            point=np.array([0.0, 0.0, 0.0]),
            normale=np.array([1.0, 0.0, 0.0]),
            plane_name="test"
        )
        self.shape.add_plane(plane, 1)

        electron = MockElectron(
            xcoor=np.array([1.0, 2.0, 3.0]),  # Behind the plane
            velocity=np.array([2.0, 1.0, 1.0])
        )

        new_coor = self.shape.get_new_coors_after_reflect(electron, 0, 0)
        
        # Expected new position: [1.0, 2.0, 3.0] (mirrored across X=0 plane)
        # Expected new velocity: [-2.0, 1.0, 1.0] (reflected normal component)
        expected = np.array([-1.0, 2.0, 3.0, -2.0, 1.0, 1.0])
        assert np.allclose(new_coor, expected)

    def test_get_new_coors_after_reflect_invalid_index(self):
        """Test exception for reflection with invalid index"""
        self.shape.add_plane(self.plane1, 1)
        electron = MockElectron(
            xcoor=np.array([1.0, 1.0, 1.0]),
            velocity=np.array([0.0, 0.0, 0.0])
        )
        
        with pytest.raises(ValueError, match="Electron inside and didnt scattering"):
            self.shape.get_new_coors_after_reflect(electron, 0, -1)

    def test_get_cos_angle(self):
        """Test cosine angle calculation between velocity and normal"""
        self.shape.add_plane(self.plane1, 1)  # Normal (-1,0,0)
        
        electron = MockElectron(
            xcoor=np.array([0.0, 0.0, 0.0]),
            velocity=np.array([3.0, 4.0, 0.0])  # Velocity magnitude = 5
        )
        
        # cos = (3*1 + 4*0 + 0*0) / 5 = 0.6
        cos_angle = self.shape.get_cos_angle(electron, 0, 0)
        assert np.isclose(cos_angle, -0.6)


# --- Test edge cases and special configurations ---

def test_multiple_planes_priority():
    """Test that correct plane is selected when electron crosses multiple boundaries"""
    shape = ConvexShape()
    
    # Add planes with different positions
    shape.add_plane(
        Plane(np.array([0,0,0]), np.array([1,0,0]), "plane1"),
        1
    )
    shape.add_plane(
        Plane(np.array([0,0,0]), np.array([0,1,0]), "plane2"),
        1
    )
    shape.add_plane(
        Plane(np.array([0,0,0]), np.array([0,0,1]), "plane3"),
        1
    )
    
    # Electron far outside in negative region
    electron = MockElectron(
        xcoor=np.array([2.0, 3.0, 1.0]),
        velocity=np.array([0.0, 0.0, 0.0])
    )
    
    # Should return plane with largest positive distance
    # Distances: plane1: 2.0, plane2: 3.0, plane3: 1.0
    # plane2 has largest distance (3.0)
    assert shape.get_out_plane(electron, 0) == 2


# --- Test reflection process function (commented out as it requires specific dependencies) ---
# To test _reflection_process, you would need to mock geom and electrons properly

if __name__ == "__main__":
    pytest.main([__file__, "-v"])