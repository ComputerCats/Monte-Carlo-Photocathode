import numpy as np
import scipy.integrate as integrate

STATUS = {'Exit': 'Exit', 'Outside': 'Outside','Inside': 'Inside','Reflect': 'Reflect', 'Died': 'Died'}

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

    def get_name(self):

        return 'HalfspaceGeom'

    def get_params(self):

        return f'point = {self.point}, normale = {self.normale}'

    def get_distance(self, other_point):

        distance = np.dot(other_point - self.point, self.normale)

        return distance

    def _is_outside(self, coor):

        if self.get_distance(coor) > 0:

            return STATUS['Outside']

        else: 
            
            return STATUS['Inside']

    def get_new_point_after_reflect(self, prev_dir, curr_point):

        result = curr_point + 2*np.array([0, 0, (self.point[2] - curr_point[2])])

        return result

    def get_new_dir_after_reflect(self, prev_dir, curr_point):

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

    def __init__(self, point_substrate, normale_substrate, point_out, normale_out, R = 1):

        self.point_substrate = point_substrate
        self.normale_substrate = normale_substrate
        self.point_out = point_out
        self.normale_out = normale_out

        self.R = R

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

    def get_new_point_after_reflect(self, direction, curr_point):

        if self.get_distance(curr_point, self.point_substrate, self.normale_substrate) > 0:

            result = curr_point[2] - 2*np.array([0, 0, (curr_point[2] - self.point_substrate[2])])

        if self.get_distance(curr_point, self.point_out, self.normale_out) > 0:

            result = curr_point + 2*np.array([0, 0, (self.point_out[2] - curr_point[2])])

        return result

    def get_new_dir_after_reflect(self, prev_dir, curr_point):

        prev_dir[1] = np.pi - prev_dir[1]

        result = prev_dir

        return result

    def _get_substrat_exit_status(self):

        substrat_reflect_or_died_status = np.random.choice([STATUS['Died'], STATUS['Reflect']], p=[1-self.R, self.R])

        return substrat_reflect_or_died_status

    def get_status(self, point):

        if self._is_outside(point) == STATUS['Outside']:

            if self.get_distance(point, self.point_out, self.normale_out) > 0:
                
                return STATUS['Exit']

            else:
                
                return self._get_substrat_exit_status()

        else:

            return STATUS['Inside']

    def get_cos_angle(self, electron): #return cos for external normal

        angle = electron.get_dir()[1]
        result = -np.cos(angle)

        return result

class Rectangular:

    #       point 2
    #/--------*       |
    #|        |       |
    #|        |       |
    #|        |       |
    #|        |       |
    #|        |       V  Z
    #*--------/
    #point 1

    #-------> X

    def __init__(self, point1, point2, substrate_normales_indx = [-1]):

        self.normales = [np.array([1, 0, 0]), np.array([-1, 0, 0]), np.array([0, 0, 1]), np.array([0, 0, -1])]
        self.point1 = point1
        self.point2 = point2
        self.substrate_normales_indx = substrate_normales_indx

    def get_name(self):

        return 'Rectangular'

    def __str__(self) -> str:
        return self.get_params()

    def get_params(self):

        return f'point1 = {self.point1}, point2 = {self.point2}, substrate_normales_indx = {self.substrate_normales_indx}'

    @staticmethod
    def get_distance(other_point, plane_point, normale):

        distance = np.dot(other_point - plane_point, normale)

        return distance

    def _is_outside(self, other_point):

        if np.dot(other_point - self.point1, self.normales[1]) > 0:

            if 1 in self.substrate_normales_indx:

                return STATUS['Died']

            else:

                return STATUS['Exit']

        if np.dot(other_point - self.point1, self.normales[2]) > 0:

            if 2 in self.substrate_normales_indx:

                return STATUS['Died']

            else:

                return STATUS['Exit']

        if np.dot(other_point - self.point2, self.normales[0]) > 0:

            if 0 in self.substrate_normales_indx:

                return STATUS['Died']

            else:

                return STATUS['Exit']

        if np.dot(other_point - self.point2, self.normales[3]) > 0:

            if 3 in self.substrate_normales_indx:

                return STATUS['Died']

            else:

                return STATUS['Exit']

        return STATUS['Inside']

    def get_outer_way(self, curr_point):

        if curr_point[0] < self.point1[0]:

            return 0

        if curr_point[0] > self.point2[0]:

            return 1

        if curr_point[2] < self.point2[2]:

            return 2

        if curr_point[2] > self.point1[2]:

            return 3

        raise ValueError('Coudnt find outer plane')

    def get_new_point_after_reflect(self, dir, curr_point):

        outer_normale_indx = self.get_outer_way(curr_point)
        if outer_normale_indx == 0:

            result = curr_point + 2*np.array([self.point1[0] - curr_point[0], 0, 0])

        if outer_normale_indx == 1:

            result = curr_point - 2*np.array([self.point2[0] - curr_point[0], 0, 0])

        if outer_normale_indx == 2:

            result = curr_point + 2*np.array([0, 0, curr_point[1] - self.point2[1]])

        if outer_normale_indx == 3:

            result = curr_point - 2*np.array([0, 0, self.point1[1] - curr_point[1]])

        return result

    def get_new_dir_after_reflect(self, prev_dir, curr_point):

        outer_normale_indx = self.get_outer_way(curr_point)

        if outer_normale_indx == 0 or outer_normale_indx == 1:

            prev_dir[0] = prev_dir[0] - np.pi

        if outer_normale_indx == 2 or outer_normale_indx == 3:

            prev_dir[1] = np.pi - prev_dir[1]

        result = prev_dir

        return result

    def get_status(self, point):

        return self._is_outside(point)

    def get_cos_angle(self, electron): #return cos for external normal

        outer_normale_indx = self.get_outer_way(electron.get_coor())

        curr_dir_phi_pheta = electron.get_dir()
        

        cos_theta_mass = np.cos(curr_dir_phi_pheta[1])
        sin_theta_mass = np.sin(curr_dir_phi_pheta[1])
        cos_psi_mass = np.cos(curr_dir_phi_pheta[0])

        vx = sin_theta_mass*cos_psi_mass
        vz = cos_theta_mass

        if outer_normale_indx == 0:

            return vx

        if outer_normale_indx == 1:

            return -vx

        if outer_normale_indx == 2:

            return vz

        if outer_normale_indx == 3:

            return -vz


class OneD:

    #normales must be vec to outside

    def __init__(self, rect1, rect2, R = 1):

        self.rect1 = rect1
        self.rect1 = rect2

        self.R = R

    def get_name(self):

        return 'OneD'

    def get_params(self):

        return f'pos_11_rect = {self.pos_11_rect}, pos_12_rect = {self.pos_12_rect}, pos_21_rect = {self.pos_21_rect}, pos_22_rect = {self.pos_22_rect}'

    @staticmethod
    def get_distance(other_point, plane_point, normale, indx_box):

        if indx_box == 1:

            distance = np.dot(other_point - plane_point, normale)

        return distance

    def _is_outside(self, coor):

        if self.get_distance(coor, self.point_substrate, self.normale_substrate) > 0 or self.get_distance(coor, self.point_out, self.normale_out) > 0:

            return STATUS['Outside']

        else: 
            
            return STATUS['Inside']

    def get_new_point_after_reflect(self, direction, curr_point):

        if self.get_distance(curr_point, self.point_substrate, self.normale_substrate) > 0:

            result = curr_point - 2*np.array([0, 0, (curr_point[2] - self.point_substrate[2])])

        if self.get_distance(curr_point, self.point_out, self.normale_out) > 0:

            result = curr_point + 2*np.array([0, 0, (self.point_out[2] - curr_point[2])])

        return result

    def get_new_dir_after_reflect(self, prev_dir, curr_point):

        prev_dir[1] = np.pi - prev_dir[1]

        result = prev_dir

        return result

    def _get_substrat_exit_status(self):

        substrat_reflect_or_died_status = np.random.choice([STATUS['Died'], STATUS['Reflect']], p=[1-self.R, self.R])

        return substrat_reflect_or_died_status

    def get_status(self, point):

        if self._is_outside(point) == STATUS['Outside']:

            if self.get_distance(point, self.point_out, self.normale_out) > 0:
                
                return STATUS['Exit']

            else:
                
                return self._get_substrat_exit_status()

        else:

            return STATUS['Inside']

    def get_cos_angle(self, electron): #return cos for external normal

        angle = electron.get_dir()[1]
        result = -np.cos(angle)

        return result