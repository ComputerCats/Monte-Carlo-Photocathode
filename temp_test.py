import numpy as np
import electron


class Test:
	def __init__(self):
	    self.N = 2
	    self.electrons = electron.Electrons(self.N)

	    test_coor = np.array([
	        [0.0, 0.0, 0.0, 1e6, 0.0, 0.0, 0], 
	        [0.0, 0.0, 0.0, 0.0, 2e6, 0.0, 0]  
	    ])
	    self.electrons.coor = test_coor.copy()
	    self.electrons.set_electron_properties(0.5)

	def test_add_energy_all(self):

	    initial_energy = self.electrons.get_E(0)
	    delta_E = 1e-18
	    print(self.electrons.get_E(0))
	    self.electrons.add_energy(np.array([delta_E, delta_E]))
	    
	    final_energy = self.electrons.get_E(0)

	    return {'res1': final_energy, 'res2': initial_energy + delta_E}


test = Test()
res = test.test_add_energy_all()
print(res)