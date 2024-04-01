import numpy as np
import scipy.integrate as integrate

STATUS = {'Exit': 'Exit', 'Outside': 'Outside','Inside': 'Inside','Reflect': 'Reflect', 'Died': 'Died'}

def trans_sphere_to_dec_norm(psi, theta):

    return np.array([np.cos(psi)*np.sin(theta), np.sin(psi)*np.sin(theta), np.cos(theta)])

def L_2_norm(func, args, min_bound, max_bound):

    norm = integrate.quad(func, min_bound, max_bound, args = args)[0]

    return norm

def get_h(point, point_electron, out_normale):

    result = np.dot(point_electron - point, out_normale)

    return result

class Plane:

    def __init__(self, point, normale) -> None:

        self.point = point
        self.normale = normale

    def get_normale(self):

        return self.normale

    def get_point(self):

        return self.point

    def __str__(self) -> str:
        
        return f'point = {self.get_point()}, normale = {self.get_normale()}'

class ConvexShape:

    def __init__(self) -> None:
        
        self.planes = []
        self.work_surfaces = []

    def add_plane(self, plane, surface_status)-> None:

        self.planes.append(plane)

        self.work_surfaces.append(surface_status)

    def get_name(self) -> str:

        return 'ConvexShape'

    def get_params(self) -> str:

        return f'planes = {self.planes}'

    def _find_exit_params_plane(self, single_electron) -> dict:
        
        curr_point = single_electron.get_coor()

        h_prev = 1000

        exit_indx = 0

        for indx, plane in enumerate(self.planes):

            h = get_h(plane.get_point(), curr_point, plane.get_normale())

            if h > 0 and h_prev > h:

                normale = plane.get_normale()
                h_prev = h
                point = plane.get_point()
                exit_indx = indx

        return {'normale': normale, 'h': h_prev, 'point': point, 'indx': exit_indx}

    def get_new_coors_after_reflect(self, single_electron) -> dict:

        curr_point = single_electron.get_coor()
        prev_dir = single_electron.get_veloicity_vector()

        reflect_params = self._find_exit_params_plane(single_electron)

        new_dir = prev_dir - 2*np.dot(reflect_params['normale'], prev_dir)*reflect_params['normale']
        new_coor = curr_point - 2*np.dot(reflect_params['normale'], curr_point - reflect_params['point'])*reflect_params['normale']

        return {'new_coor': new_coor, 'new_dir': new_dir}
        
    def get_status(self, single_electron) -> str:

        curr_point = single_electron.get_coor()

        for plane in self.planes:

            if get_h(plane.get_point(), curr_point, plane.get_normale()) > 0:

                exit_params = self._find_exit_params_plane(single_electron)

                return self.work_surfaces[exit_params['indx']]

        return STATUS['Inside']

    def get_cos_angle(self, electron): #return cos for external normal

        reflect_params = self._find_exit_params_plane(electron)

        result = np.dot(electron.get_veloicity_vector(), reflect_params['normale'])/electron.get_veloicity()

        return result
