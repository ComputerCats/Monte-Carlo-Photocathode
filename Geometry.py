import numpy as np
import scipy.integrate as integrate

STATUS = {'Exit': 'Exit', 'Outside': 'Outside','Inside': 'Inside','Reflect': 'Reflect','~reflect': '~reflect', 'Died': 'Died'}

def trans_sphere_to_dec_norm(psi, theta):

    return np.array([np.cos(psi)*np.sin(theta), np.sin(psi)*np.sin(theta), np.cos(theta)])

def L_2_norm(func, args, min_bound, max_bound):

    norm = integrate.quad(func, min_bound, max_bound, args = args)[0]

    return norm

class HalfspaceGeom:

    #normales must be vec to outside

    def __init__(self, point, normale):

        self.point = point
        self.normale = normale

        self.reflection_coef = 1

    def get_name(self):

        return 'HalfspaceGeom'

    def get_params(self):

        return f'point = {self.point}, normale = {self.normale}'

    def set_absorption(self, reflection_coef):

        self.reflection_coef = reflection_coef

    def get_distance(self, other_point):

        distance = np.dot(other_point - self.point, self.normale)

        return distance

    def _is_outside(self, coor):

        if self.get_distance(coor) > 0:

            return STATUS['Outside']

        else: 
            
            return STATUS['Inside']

    def get_new_point_after_reflect(self, curr_point, direction):

        result = curr_point + 2*np.array([0, 0, (self.point[2] - curr_point[2])])

        return result

    def get_new_dir_after_reflect(self, prev_dir):

        prev_dir[1] = np.pi - prev_dir[1]

        result = prev_dir

        return result

    def get_status(self, point):

        if self._is_outside(point) == STATUS['Outside']:
            
            return STATUS['Exit']

        else:

            return STATUS['Inside']

    def get_cos_angle(self, electron): #return cos for external normal

        angle = electron.get_dir()[1]
        result = -np.cos(angle)

        return result

class PlateGeom:

    #normales must be vec to outside

    def __init__(self, point_substrate, normale_substrate, point_out, normale_out):

        self.point_substrate = point_substrate
        self.normale_substrate = normale_substrate
        self.point_out = point_out
        self.normale_out = normale_out

    def get_name(self):

        return 'PlateGeom'

    def get_params(self):

        return f'point_substrate = {self.point_substrate}, normale_substrate = {self.normale_substrate}, point_out = {self.point_out}, normale_out = {self.normale_out}'

    @staticmethod
    def get_distance(other_point, plane_point, normale):

        distance = np.dot(other_point - plane_point, normale)

        return distance


    def _is_outside(self, coor):

        if self.get_distance(coor, self.point_substrate, self.normale_substrate) > 0 or self.get_distance(coor, self.point_out, self.normale_out) > 0:

            return STATUS['Outside']

        else: 
            
            return STATUS['Inside']

    def get_new_point_after_reflect(self, curr_point, direction):

        if self.get_distance(curr_point, self.point_substrate, self.normale_substrate) > 0:

            result = curr_point - 2*np.array([0, 0, (curr_point[2] - self.point_substrate[2])])

        if self.get_distance(curr_point, self.point_out, self.normale_out) > 0:

            result = curr_point + 2*np.array([0, 0, (self.point_out[2] - curr_point[2])])

        return result

    def get_new_dir_after_reflect(self, prev_dir):

        prev_dir[1] = np.pi - prev_dir[1]

        result = prev_dir

        return result

    def get_status(self, point):

        if self._is_outside(point) == STATUS['Outside']:

            if self.get_distance(point, self.point_out, self.normale_out) > 0:
                
                return STATUS['Exit']

            else:
                
                return STATUS['Reflect']

        else:

            return STATUS['Inside']

    def get_cos_angle(self, electron): #return cos for external normal

        angle = electron.get_dir()[1]
        result = -np.cos(angle)

        return result
