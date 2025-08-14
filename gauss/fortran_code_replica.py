"""
[MJP Comment]: Add a description of the module here.
"""
import pytest
import numpy as np
from gauss import math_fucntions as mf


# [MJP Comment]: If I understand correctly what these position_vector* and unit_vector* and observation_times variables are,
# [MJP Comment]: these quantities should be moved into the `test_fortran_code_replica` module, and then you should add a
# [MJP Comment]: description (in `test_fortran_code_replica.py`) of what these variables represent and where they come from.
postion_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
postion_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
postion_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
postion_matrix = np.matrix([[postion_vector1[0],postion_vector2[0],postion_vector3[0]],[postion_vector1[1],postion_vector2[1],postion_vector3[1]],[postion_vector1[2],postion_vector2[2],postion_vector3[2]]])


unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
unit_matrix = np.matrix([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],[unit_vector1[1],unit_vector2[1],unit_vector3[1]],[unit_vector1[2],unit_vector2[2],unit_vector3[2]]])

observation_times = [58577.489970740738,58583.545550740739,58590.545550740739]

def get_tau_values(time_obs,sqrt_mu):
    """
    [MJP Comment]: Add a description of the function here (and do the same for all other functions).
    """
    tau_1 = sqrt_mu * (time_obs[0] - time_obs[1]) 
    tau_3 = sqrt_mu * (time_obs[2] - time_obs[1])
    tau = tau_3 - tau_1
    return (tau,tau_1,tau_3)

def make_a_b(tau,tau1,tau3):
    a_1 = tau3 / tau
    b_1 = a_1 * (tau**2 - tau3**2) / 6
    
    a_3 = -tau1 / tau
    b_3 = a_3 * (tau**2 - tau1**2) / 6
    return (a_1,b_1,a_3,b_3)


def get_d_and_pos_variables(inv_matrix,dot_b,dot_a,obs_matrix,unit_vect):
    d_1 = float(inv_matrix[1,0]) * float(dot_a[0,0]) + float(inv_matrix[1,1])* float(dot_a[1,0]) + float(inv_matrix[1,2]) * float(dot_a[2,0]) # a2star

    dd_2 = float(inv_matrix[1,0] * dot_b[0,0] + inv_matrix[1,1] * dot_b[1,0] + inv_matrix[1,2] * dot_b[2,0]) # b2star
    
    r2_pos = obs_matrix[0,1]**2 +obs_matrix[1,1]**2 + obs_matrix[2,1]**2   #r22, r_2 magnitude vector square

    dot_prod_r2_comp = unit_vect[0,1]*obs_matrix[0,1] + unit_vect[1,1]*obs_matrix[1,1] + unit_vect[2,1]*obs_matrix[2,1] # s22r2
    return (d_1, dd_2, float(r2_pos),float(dot_prod_r2_comp))


def generate_Cs(d_1,d_2,r2pos,dotr2):
    c_8 = 1
    c_6 = -(d_1**2)-r2pos-(2*d_1*dotr2)
    c_3 = -(2*d_2*(d_1+dotr2))
    c = -(d_2**2)
    return (c_8,c_6,c_3,c)


def calculate_rhos(a1,b1,a3,b3,r2m3,pos_matrix, inv_matrix): 
    c_1 = a1 + b1*r2m3
    c_3 = a3 + b3 *r2m3      
    matrix_c = [[c_1],[-1],[c_3]]
    gcap = mf.matrix_dot_prod(pos_matrix,matrix_c)
    crhom = mf.matrix_dot_prod(inv_matrix,gcap)
    p_1 = -(crhom[0,0]/c_1)
    p_3 = -(crhom[2,0]/c_3)
    return (float(p_1),float(p_3),c_1,c_3,gcap,crhom,r2m3)


def find_rho_and_r(root,a1,b1,a3,b3,pos_matrix,inv_matrix,unit_vect):
        r2m3 = 1/(root**3)
        rhos = calculate_rhos(a1,b1,a3,b3,r2m3,pos_matrix,inv_matrix)
        rhos = [rhos[0],rhos[1]]
        r_1 = pos_matrix[0:,0] + rhos[0]*unit_vect[:,0]
        r_3 = pos_matrix[0:,2] + rhos[-1]*unit_vect[:,2]
        r = [r_1, r_3] # matrix of 3x2
        return(rhos,r)


def gauss_method_8th(obs_postions,times_obs,obj_unit_vector):
    mu = 1.7202098949957226E-002 # mue
    inv_unit_vector = mf.take_inverse_matrix(obj_unit_vector.copy())
    tau, tau_1, tau_3 =get_tau_values(times_obs,mu)
    
    a_1, b_1, a_3, b_3 = make_a_b(tau,tau_1,tau_3)
    
    matrix_a = [[a_1],[-1],[a_3]]
    matrix_b = [[b_1],[0],[b_3]]

    pos_dot_a = mf.matrix_dot_prod(obs_postions,matrix_a)
    pos_dot_b = mf.matrix_dot_prod(obs_postions,matrix_b)
    
    d_1, dd_2, r2_pos,dot_prod_r2_comp=get_d_and_pos_variables(inv_unit_vector,pos_dot_b,pos_dot_a,obs_postions,obj_unit_vector)
    
    c_8,c_6,c_3,c = generate_Cs(d_1,dd_2,r2_pos,dot_prod_r2_comp)
    
    roots = mf.eight_equation_four_coeff(c,c_3,c_6,c_8)
    results = []
    for i in range(len(roots)):
        rho, r = find_rho_and_r(roots[i],a_1,b_1,a_3,b_3,obs_postions,inv_unit_vector,obj_unit_vector)
        results.append([roots[i],rho,r])
    return results


# print(gauss_method_8th(postion_matrix,observation_times,unit_matrix))
