import electron

N_electrons = 10
N_iterations = 10
E_g = 1.2
E_a = 0.7
gamma = 2.5
delta_E = 0.015

task = electron.Electrons(N_electrons, N_iterations, E_g, E_a, gamma, delta_E)
task.initial_phase_prostr()
task.run_simulation()






