import emtl
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append(r'C:\Users\Mikhail\YandexDisk\Kintech\Pwr\Pwr')
sys.path.append(r'C:\Users\Mikhail\YandexDisk\Kintech\Pwr\Pwr\models')

import MatLib as materials
import Pwr as pwr
import Halfspace as halfspace

time = 100
second_time = 60

way_to = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\fdtd\halfspace\results'
way_from = r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\Monte-Carlo-Photocathode\fdtd\halfspace'

task = pwr.Pwr('test_halfspace', way_to, 'test')

model = halfspace.Halfspace()
model.chose_blockId(0)
model.set_material(materials.get_CsSbK2(), 'K2CsSb')
model.set_dx(0.003)
model.set_geom_params(zsize = 1.024)
model.set_signal(theta = 0)
model.set_PML_type('NPML')
#model.set_semi(materials.get_CsSbK2(), 0.024)
#model.set_photo_semi_detector('semi_detector', 0.003)
#model.set_detectors_reflection()
model.set_fourie()

print(model.get_model_params())

task.chose_model(model)
task.run_experiment(time, second_time)
task.save_result(way_from)