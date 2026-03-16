import numpy as np

STATUS = {'Exit': 1, 'Died': 0, 'Reflect': -1}

def _get_indx_min_positive_value(hs):

    indx_min_max_pos_value = np.argmax(hs)
    min_max_pos_value = hs[indx_min_max_pos_value]

    for indx, element in enumerate(hs):

        if element >= 0 and element <= min_max_pos_value:

            min_max_pos_value = element
            indx_min_max_pos_value = indx

    return indx_min_max_pos_value

class Plane:

    def __init__(self, point, normale, plane_name) -> None:

        self.point = point
        self.normale = normale
        self.plane_name = plane_name

    def get_normale(self):

        return self.normale

    def get_point(self):

        return self.point

    def get_distance(self, coor):

        result = np.dot(coor - self.point, self.normale)
        
        return result

    def __str__(self) -> str:
        
        return f'point = {self.get_point()}, normale = {self.get_normale()}, plane_name = {self.plane_name}'

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

        res_str = ''

        for plane in self.planes:

            res_str += f' {str(plane)}'

        return res_str

    def get_out_plane(self, electrons, el_indx: int) -> int:

        xcoor = electrons.get_xcoor(el_indx)
        hs = [plane.get_distance(xcoor) for plane in self.planes]
        max_h = np.max(hs)

        if max_h < 0:

            return -1

        else:

            return _get_indx_min_positive_value(hs)

    def get_type_exit_plane(self, indx_plane: int) -> int:

        if indx_plane == -1:

            raise ValueError('Electron inside and cannt escape')

        return self.work_surfaces[indx_plane]

    def get_new_coors_after_reflect(self, electrons, el_indx: int, indx_out: int):

        if indx_out == -1:

            raise ValueError('Electron inside and didnt scattering')

        exit_plane = self.planes[indx_out]
        normale = exit_plane.get_normale()
        point = exit_plane.get_point()

        curr_point = electrons.get_xcoor(el_indx)
        prev_dir = electrons.get_velocity(el_indx)

        new_dir = prev_dir - 2*np.dot(normale, prev_dir)*normale
        new_coor = curr_point - 2*np.dot(normale, curr_point - point)*normale

        return np.hstack((new_coor, new_dir))
        
    def get_cos_angle(self, electrons, el_indx: int, indx_out: int): #return cos for external normal

        module_velocity = electrons.get_module_velocity(el_indx)

        if module_velocity <= 0:

            raise ValueError('Module of velocity is zero')

        result = np.dot(electrons.get_velocity(el_indx), self.planes[indx_out].get_normale())/module_velocity

        return result
