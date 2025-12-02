'''
This file is a replica of the fortran code used to produce a result for the gauss method 
using the eight degree equation, the content and procedures and variable names are 
the same as the one presented within the fortran code, and the results of this 
code are able to replicate the fortran results up to the 12th decimal place. this code replica
was made to get a further understanding into the eight degree gauss method, a more 'independent' 
version that is not a replica can be found in the file gauss_method.py'''
import numpy as np
from gauss import math_fucntions as mf
from gauss import demo 

def get_tau_values(time_obs: list,rootgm: float = 1.7202098949957226E-002):
    """Provided with a set of three times of observations and mue square value
    it produces the tau,tau_1,tau_3 necessary for the gauss method

    Args:
        time_obs (list): a set of three observation times in MJD
        rootgm (float): the value of mue square: 1.7202098949957226E-002

    Returns:
        tuple: the resulting tau values as a tuple: (tau,tau_1,tau_3) 
    """
    tau_1 = rootgm * (time_obs[0] - time_obs[1]) 
    tau_3 = rootgm * (time_obs[2] - time_obs[1])
    tau13 = tau_3 - tau_1
    return (tau13,tau_1,tau_3)

def make_a_b(tau13: float,tau1: float,tau3:float):
    """Provided three tau values (tau,tau_1,tau_3), calculated with the observation time
    it produces the a_1,b_1,a_3,b_3 necessary for the gauss method.

    Args:
        tau (float): tau value calculated using the time observations mathematical 
        procedure is in function get_tau_values.
        tau1 (float): tau1 value calculated using the time observations mathematical 
        procedure is in function get_tau_values.
        tau3 (float): tau3 value calculated using the time observations mathematical 
        procedure is in function get_tau_values.

    Returns:
        tuple: The resulting values as a tuple (a_1,b_1,a_3,b_3)
    """
    a_1 = tau3 / tau13
    b_1 = a_1 * (tau13**2 - tau3**2) / 6
    a_3 = -tau1 / tau13
    b_3 = a_3 * (tau13**2 - tau1**2) / 6
    return (a_1,b_1,a_3,b_3)

def get_d_and_pos_variables(sinv0:np.matrix,rb:np.matrix,ra:np.matrix,xt,esse0:np.matrix):
    """Computes scalar parameters (a2star, b2star, r22, s2r2) used in Gauss method.
    Args:
        sinv0 (np.matrix): 3x3 inverse matrix of the unit vector matrix esse 0 
        rb (np.matrix): a matrix of 3x1 containing the dot product between xt and [b1 0 b3] 
        ra (np.matrix): a matrix of 3x1 containing the dot product between xt and [a1 1 a3]
        xt (np.matrix): 3x3 matrix representing the observer position matrix where each column represents a vector
        esse0 (np.matrix): 3x3 matrix, containing the unit vector of the observation

    Returns:
        tuple: a tuple of length four containing the resulting values (a2star, b2star, r22,s2r2)
    """
    a2star = float(sinv0[1,0]) * float(ra[0,0]) + float(sinv0[1,1])* float(ra[1,0]) + float(sinv0[1,2]) * float(ra[2,0]) 
    b2star = float(sinv0[1,0] * rb[0,0] + sinv0[1,1] * rb[1,0] + sinv0[1,2] * rb[2,0]) 
    r22 = xt[0,1]**2 +xt[1,1]**2 + xt[2,1]**2   
    s2r2 = esse0[0,1]* xt[0,1] + esse0[1,1]* xt[1,1] + esse0[2,1]*xt[2,1] 
    return (a2star, b2star, float(r22),float(s2r2))

def generate_Cs(a2star:float,b2star: float,r22: float,s2r2: float):
    """Provided with a2star:float,b2star: float,r22: float,s2r2 calculated from function get_d_and_pos_variables
    calculates the coefficient values for the eight degree equation in the gauss method.

    Args:
        a2star (float): a float value calculated from function get_d_and_pos_variables
        b2star (float): a float value calculated from function get_d_and_pos_variables
        r22 (float): magnitude of the observer position matrix xt
        s2r2 (float): a float value calculated from function get_d_and_pos_variables

    Returns:
        tuple: a tuple containing four coefficient that will be used for the eight degree method (c_8,c_6,c_3,c)
    """
    c_8 = 1
    c_6 = -(a2star**2)-r22-(2*a2star*s2r2)
    c_3 = -(2*b2star*(a2star + s2r2))
    c = -(b2star**2)
    return (c_8,c_6,c_3,c)

def calculate_rhos(a1: float,b1: float,a3: float,b3: float,r2m3: float,xt:np.matrix, sinv0:np.matrix):
    """This function calculates the rhos values given, a1,b1,a3,b3,r2m3,xt, sinv0

    Args:
        a1 (float): a1 value calculated from function make_a_b
        b1 (float): b1 value calculated from function make_a_b
        a3 (float): a3 value calculated from function make_a_b
        b3 (float): b3 value calculated from function make_a_b
        r2m3 (float): Value calculated from 1/(root**3)
        xt (np.matrix): the 3x3 observer position matrix
        sinv0 (np.matrix): the 3x3 inverse of the unit vector matrix esse0

    Returns:
        tuple: where the first two provided values are rho1 and rho3 and the other provides values are for testing
        returned values are in the following order rho_1,rho_3,c_1,c_3,gcap,crhom,r2m3
    """
    c_1 = a1 + b1*r2m3
    c_3 = a3 + b3 *r2m3      
    temp_matrix = [[c_1],[-1],[c_3]]
    gcap = mf.matrix_dot_prod(xt,temp_matrix)
    crhom = mf.matrix_dot_prod(sinv0,gcap)
    rho_1 = -(crhom[0,0]/c_1)
    rho_3 = -(crhom[2,0]/c_3)
    return (float(rho_1),float(rho_3),c_1,c_3,gcap,crhom,r2m3)

def find_rho_and_r(root: float,a1: float,b1: float,a3: float,b3: float,xt:np.matrix,sinv0:np.matrix,esse0: np.matrix):
    """This function depends on calculate_rhos. 
    This function provided the root, a1, b1, a3, b3, observer position matrix xt, unit vector esse0, and the inverse of the 
    unit vector sinv0, to calculate the rho values, and the the position of the asteroid xp 
    
    Args:
        root (float): The current root in question solved from the eight degree equation
        a1 (float): a1 value calculated from function make_a_b
        b1 (float): b1 value calculated from function make_a_b
        a3 (float): a3 value calculated from function make_a_b
        b3 (float): b3 value calculated from function make_a_b
        xt (np.matrix): 3x3 observer position matrix
        sinv0 (np.matrix): a  3x3 inverse of the unit vector matrix esse0
        esse0 (np.matrix): the 3x3 unit vector matrix
    """
    r2m3 = 1/(root**3)
    rhos = calculate_rhos(a1,b1,a3,b3,r2m3,xt,sinv0)
    rhos = [rhos[0],rhos[1]]
    xp_1 = xt[0:,0] + rhos[0]*esse0[:,0]
    xp_3 = xt[0:,2] + rhos[-1]*esse0[:,2]
    xp = [xp_1, xp_3] # matrix of 3x2
    return(rhos,xp)

def gauss_method_8th(xt,times_obs:list,esse0):
    """main function that call on and depends on functions matrix_dot_prod, 
    take_inverse_matrix, eight_equation_four_coeff from the math_fucntions.py 
    file, and on local functions get_tau_values, make_a_b, 
    get_d_and_pos_variables, generate_Cs and find_rho_and_r. 
    
    This function provided the observer position xt, tree times of observations in time_obs, and the unit vector of observation as esse0, it calculates the gauss method that 
    produces multiple roots, and for each of the roots calculates the position of the asteroid xp, and the rho values

    Args:
        xt (list or np.matrix): _description_
        times_obs (list): a list of length 3 containing 3 set of times when the observation happen in MJD
        esse0 (list or np.matrix): 3x3 Unit vector matrix

    Returns:
        list: a list of lengths of 3 or less where each element is a sublist containing 3 elements, 
        the root calculated in the eight degree equation, the the rhos 
        values and xp the position.
    """
    xt = np.matrix(xt)
    esse0 = np.matrix(esse0)
    mu = 1.7202098949957226E-002 # mu
    sinv0 = mf.take_inverse_matrix(esse0.copy())
    tau13, tau_1, tau_3 =get_tau_values(times_obs,mu)
    
    a_1, b_1, a_3, b_3 = make_a_b(tau13,tau_1,tau_3)
    
    matrix_a = [[a_1],[-1],[a_3]]
    matrix_b = [[b_1],[0],[b_3]]
    ra = mf.matrix_dot_prod(xt,matrix_a)
    rb = mf.matrix_dot_prod(xt,matrix_b)
    a2star, b2star, r22,s2r2p = get_d_and_pos_variables(sinv0,rb,ra,xt,esse0)
    c_8,c_6,c_3,c = generate_Cs(a2star, b2star, r22,s2r2p)
    roots = mf.eight_equation_four_coeff(c,c_3,c_6,c_8)
    results = []
    for i in range(len(roots)):
        rho, xp = find_rho_and_r(roots[i],a_1,b_1,a_3,b_3,xt,sinv0,esse0)
        results.append([roots[i],rho,xp])
    return results

if __name__ == "__main__":
    demo_result = gauss_method_8th(demo.xt,demo.observation_times,demo.esse0)
    print("this are the roots produced by the gauss method")
    for i in range(len(demo_result)):
        root = demo_result[i][0]
        print("root: ", root)
