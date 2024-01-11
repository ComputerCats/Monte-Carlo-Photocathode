import numpy as np

#semiconductor propeties
E_G_CONST = 1.2
T_SEMICONDUCTOR = 300
DELTA_EPSILON = (1/5-1/8.2)#(1/6-1/14)
DELTA_E_POP_PLUS = 0.027
DELTA_E_POP_MINUS = 0.027
EFFECTIVE_MASS = 0.12

#HELP FOR TAU POP
CONST_PART_TAU_POP = (1.6)**3*DELTA_E_POP_PLUS/(4*np.pi*8.85)

def _Bose_Distr(E_phonon):

    return 1/(np.exp(1.6*1e4/(1.38*T_SEMICONDUCTOR)*E_phonon) - 1)


def tau_POP_plus(E):

    logarithm_part = np.log((1+np.sqrt(1+DELTA_E_POP_PLUS/E))/np.abs(1-np.sqrt(1+DELTA_E_POP_PLUS/E)))
    veloisity_reaction = 100*CONST_PART_TAU_POP*_Bose_Distr(DELTA_E_POP_PLUS)*DELTA_EPSILON*np.sqrt(EFFECTIVE_MASS/(2*E*1.6))*logarithm_part

    return 1/veloisity_reaction

def tau_POP_minus(E):

    logarithm_part = np.log((1+np.sqrt(1-DELTA_E_POP_PLUS/E))/np.abs(1-np.sqrt(1-DELTA_E_POP_PLUS/E)))
    veloisity_reaction = 100*CONST_PART_TAU_POP*(_Bose_Distr(DELTA_E_POP_PLUS)+1)*DELTA_EPSILON*np.sqrt(EFFECTIVE_MASS/(2*E*1.6))*logarithm_part

    return 1/veloisity_reaction