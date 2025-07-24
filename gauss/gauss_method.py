import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from gauss import math_fucntions as mf



mue = 1.7202098949957226E-002 **2 # Only global variable

def calculate_delta_t(times):
    d_t1 =times[0] - times[1] 
    d_t3 = times[2] - times[1]
    d_t = d_t3 - d_t1
    return (d_t1,d_t3,d_t)

def calculate_a_and_b(d_t1,d_t3,d_t):
    a_1 = d_t3 / d_t
    b_1 = (d_t3 * (d_t**2 - d_t3**2)) / (6 * d_t)
    a_3 = - d_t1 / d_t
    b_3 = - (d_t1 * (d_t**2 - d_t1 **2)) / (6 * d_t)
    return (a_1,b_1,a_3,b_3)

def calculate_values_for_cs(a_1,b_1,a_3,b_3,position_R,rho,B_matrix):
    d_1 = float(B_matrix[1,0]) * a_1 - float(B_matrix[1,1]) + float(B_matrix[1,2]) * a_3
    d_2 = float(B_matrix[1,0]) * b_1 + float(B_matrix[1,2]) * b_3
    magnitude_R_2 = np.sqrt((float(R[0,1])** 2 + float(R[1,1])** 2 + float(R[2,1])** 2))
    temp_matrix1 = [[rho[0,1],rho[1,1],rho[2,1]]] # vector components of p_2
    temp_matrix2 = [[position_R[0,1]],[position_R[1,1]],[position_R[2,1]]] # vecotr componets of R_2
    pos_2_dot_p2 = mf.matrix_dot_prod(temp_matrix1,temp_matrix2) # dot product between P_2 * R_2
    pos_2_dot_p2 = float(pos_2_dot_p2[0,0])
    return (d_1,d_2,float(magnitude_R_2),pos_2_dot_p2)

def calculate_cs(a_1,b_1,a_3,b_3,position_R,rho,B_matrix):
    d_1,d_2,magnitude_R_2,pos_2_dot_p2= calculate_values_for_cs(a_1,b_1,a_3,b_3,position_R,rho,B_matrix)
    c = float(d_2 ** 2 * mue**2)
    c_3 =float(mue * 2 * (d_1 * d_2 + d_2 * pos_2_dot_p2))
    c_6 = float(d_1**2 +magnitude_R_2**2 + 2 * d_1 *pos_2_dot_p2)
    c_8 = -1    # it can be -1 and, or multiply all the other coefficients instead of this one negative
    return(c,c_3,c_6,c_8)

def calculate_rhos_and_r(sigma,B_matrix,a_1,b_1,a_3,b_3,R):
    mag_p1 = (float(B_matrix[0,0]) * a_1 - float(B_matrix[0,1]) + float(B_matrix[0,2]) * a_3 + (float(B_matrix[0,0]) * b_1 + float(B_matrix[0,2]) * b_3) * sigma)/((-1)*(a_1+b_1*sigma))
    mag_p2 = (float(B_matrix[1,0]) * a_1 - float(B_matrix[1,1]) + float(B_matrix[1,2]) * a_3 + (float(B_matrix[1,0]) * b_1 + float(B_matrix[1,2]) * b_3) * sigma)
    mag_p3 = (float(B_matrix[2,0]) * a_1 - float(B_matrix[2,1]) + float(B_matrix[2,2]) * a_3 + (float(B_matrix[2,0]) * b_1 + float(B_matrix[2,2]) * b_3) * sigma)/ ((-1)* (a_3+b_3*sigma))
    r1 = R[:,0] + mag_p1*p[:,0]
    r2 = R[:,1] + mag_p2*p[:,1]
    r3 = R[:,2] + mag_p3*p[:,2]
    return[['placeholder'], [mag_p1,mag_p2,mag_p3],[r1,r2,r3]]


postion_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
postion_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
postion_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
postion_matrix = np.matrix([[postion_vector1[0],postion_vector2[0],postion_vector3[0]],[postion_vector1[1],postion_vector2[1],postion_vector3[1]],[postion_vector1[2],postion_vector2[2],postion_vector3[2]]])


unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
unit_matrix = np.matrix([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],[unit_vector1[1],unit_vector2[1],unit_vector3[1]],[unit_vector1[2],unit_vector2[2],unit_vector3[2]]])

# we know 3x3 matrices r(persons position) and P(unit vector) ans also t(time) list of lenght 3

def gauss_method(observation_times,p,R):
    # matric B is inverse of unit vector p
    inverse_p = mf.take_inverse_matrix(p)
    # then taking the inverse matrix and multiplying it by R
    B_matrix = mf.matrix_dot_prod(inverse_p,R)
    delta_t1,delta_t3,delta_t= calculate_delta_t(observation_times)
    a1,b1,a3,b3 = calculate_a_and_b(delta_t1,delta_t3,delta_t)
    c,c3,c6,c8= calculate_cs(a1,b1,a3,b3,R,p,B_matrix)
    roots = mf.eight_equation_four_coeff(c,c3,c6,c8)
    results = []
    for r in roots:
        sigma = mue/(r)**3
        info_root = calculate_rhos_and_r(sigma,B_matrix,a1,b1,a3,b3,R)
        info_root[0] = r
        results.append(info_root)
    return results

p = unit_matrix.copy()
R = postion_matrix.copy()
observation_times = [58577.489970740738,58583.545550740739,58590.545550740739]
# gauss_method(observation_times,p,R)