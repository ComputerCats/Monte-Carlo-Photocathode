import numpy as np
#import ivutils
#import emtl
import scipy.integrate as integrate
#a geometry class for Monte Carlo simulation. a convex body is defined by the intersection of planes 

STATUS = {'Exit': 'Exit','Inside': 'Inside','Reflect': 'Reflect','~reflect': '~reflect'}

def trans_sphere_to_dec_norm(psi, theta):

    return np.array([np.cos(psi)*np.sin(theta), np.sin(psi)*np.sin(theta), np.cos(theta)])

def L_2_norm(func, args, min_bound, max_bound):

    norm = integrate.quad(func, min_bound, max_bound, args = args)[0]

    return norm

'''
def numpy_vector_to_emtl(vector):

    result = emtl.Vector_3(vector[0], vector[1], vector[2])

    return result

class BasicGeom:
    
    def __init__(self):

        self.planes = ivutils.Plane3Vector()
        
    def add_plane(self, plane):

        self.planes.push_back(plane)

    def init_polygon(self):

        self.polygon = emtl.Polyhedron_3(self.planes)

    def is_in(self, point):

        return self.polygon.TestPoint(point)

    def get_projection(self, point):

        result = self.polygon.SurfProject(point)

        return result
'''
class HalfspaceGeom:

    def __init__(self, point, normale):

        self.point = point
        self.normale = normale

    def get_distance(self, other_point):

        distance = np.dot(other_point - self.point, self.normale)

        return distance

    def is_exit(self, other_point):

        if other_point[2] < self.point[2]:

            return 'Exit'

        else: 

            return 'Inside'

    def is_reflect(self, point):

        if not True:
            return 'Reflect'
        else:
            return '~reflect'

    #make_reflation

    def reflect(self, electron_veloisity):

        electron_veloisity[1] = np.pi - electron_veloisity[1]

        return electron_veloisity

    def get_status(self, point):

        if self.is_reflect(point) == '~reflect':

            return self.is_exit(point)

        else: 

            return 'Reflect'

    def get_cos_angle(self, electron): #return cos for external normal
        angle = electron[4]
        result = -np.cos(angle)

        return result
