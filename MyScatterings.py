import numpy as np
import pandas as pd
from scipy import interpolate

#semiconductor propeties
E_G_CONST = 1.2
T_SEMICONDUCTOR = 300
DELTA_EPSILON = (1/5-1/8.2)#(1/6-1/14)
DELTA_E_POP_PLUS = 0.027
DELTA_E_POP_MINUS = 0.027
EFFECTIVE_MASS = 0.23

#HELP FOR TAU POP
CONST_PART_TAU_POP = (1.6)**3*DELTA_E_POP_PLUS/(4*np.pi*8.85)

L_E_E = pd.read_csv(r'C:\Users\Mikhail\YandexDisk\Kintech\Projects\Science\monte-catlo\mat_prop\data\K2CsSb\l_e_e_K2CsSb.csv', sep = '; ', header = None)
L_E_E.iloc[:, 1] = L_E_E.iloc[:, 1]/10000  #in mu
FUNC_l_e_e = interpolate.interp1d(L_E_E.iloc[:, 0], L_E_E.iloc[:, 1])

def l_e_e(E):

    if E < E_G_CONST:

        return 1000

    if E - E_G_CONST < 0.09:

        return 0.05

    if E - E_G_CONST > 1.9:

        return 0.012
        
    return FUNC_l_e_e(E - E_G_CONST)
    
def _Bose_Distr(E_phonon):

    return 1/(np.exp(1.6*1e4/(1.38*T_SEMICONDUCTOR)*E_phonon) - 1)


def tau_POP_plus(E):
    '''
    logarithm_part = np.log((1+np.sqrt(1+DELTA_E_POP_PLUS/E))/np.abs(1-np.sqrt(1+DELTA_E_POP_PLUS/E)))
    veloisity_reaction = 100*CONST_PART_TAU_POP*_Bose_Distr(DELTA_E_POP_PLUS)*DELTA_EPSILON*np.sqrt(EFFECTIVE_MASS/(2*E*1.6))*logarithm_part
    
    return 1/veloisity_reaction
    '''
    return 150
    

def tau_POP_minus(E):
    '''
    logarithm_part = np.log((1+np.sqrt(1-DELTA_E_POP_PLUS/E))/np.abs(1-np.sqrt(1-DELTA_E_POP_PLUS/E)))
    veloisity_reaction = 100*CONST_PART_TAU_POP*(_Bose_Distr(DELTA_E_POP_PLUS)+1)*DELTA_EPSILON*np.sqrt(EFFECTIVE_MASS/(2*E*1.6))*logarithm_part

    return 1/veloisity_reaction
    '''
    return 20
    

class Semiconductor:

    def __init__(self, E_a, E_g, effective_mass):

        self.E_a = E_a
        self.E_g = E_g
        self.effective_mass = effective_mass

    def get_E_a(self):

        return self.E_a

    def get_E_g(self):

        return self.E_g

    def get_effective_mass(self):

        return self.effective_mass
