"""This file contains the functions used to perform the eighth-degree equation  
of the Gauss method. It also includes functions to calculate the orbital  
elements based on the resulting roots. (File still in progress.)  
"""
from gauss import math_fucntions as mf
from astropy.time import Time
import numpy as np
import sys
import pyorb
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
module_path = "/home/mila/mpc-software/wis/wis" # contains path for wis module only usuable on linux
if module_path not in sys.path:
    sys.path.append(module_path)
np.set_printoptions(threshold=np.inf,precision=25) # used to visualize all values in a numpy matrix

mue = 1.7202098949957226E-002 **2 # global variable

def flat_variable(variable):
    """Given a variable (which may be a single number or a list containing multiple elements), 
    this function ensures the variable is placed inside a list and then flattens it into a 
    one-dimensional structure. The process is stable, meaning the original order of elements 
    is preserved, but any nested lists or brackets are fully flattened into a single long list.

    Args:
        variable (list or int/float): The variable to be flattened.

    Returns:
        numpy.ndarray: A one-dimensional NumPy array containing all elements from the input variable.
    """
    try: 
        variable = np.atleast_1d([variable]) # Result keep [] brackets and changed the shape of it
        variable = variable.flatten() # takes away the brackets
    except:
        print("Error, on function flat_variable, file gauss_method.py")
        print("make sure the provided values can be converted to numpy arrays and be flatten.")
        print("Provided variable: ",variable)
        raise ValueError
    return variable

def mjd_to_jd(date):
    """converts a modified julian date to, Standard julian date
    Args:
        date (float or list): The date to be converted can be as float 
        or a list where all the elements represent a list

    Returns:
        date: the date in Julian date format
    """
    date = flat_variable(date)
    date1 = date + 2400000.5
    return np.array(date1)

def position_obs_wis(mjd,obscode): # uses Wis.py
    """Uses Wis library to calculate the position of the observatory in a Heliocentric equatorial frame. 
    Provided the mjd and observatory code. the mjd can be provided as sets of dates, however the 
    observatory code can only be a single string

    Args:
        mjd (float or list): list of dates of the observation 
        obscode (string): Observatory code of the observation

    Returns:
        np.matrix: a 3x3 matrix per each of the given times 
        (such as provided a list of 4 times, it returns a list with 4, 3x3 matrices corresponding to each time)
    """
    # flattening the parameters 
    mjd = flat_variable(mjd)
    jd_times = mjd_to_jd(mjd)
    
    # changing the type of the provided dates into times
    times = Time(jd_times, format='jd', scale='tdb',)
    
    # importing wis to calculate the position of the observer
    import wis # imported withing the function, because  importing takes some runtime
    W = wis.wis(str(obscode), times)
    observer_posns = W.hXYZ
    return np.matrix(observer_posns)

def unit_vector_from_ra_dec(ra,dec):
    """Provided the declination and right ascension calculates the unit vector

    Args:
        ra (float or list): Right Ascension as a list of multiple floats or as a single float
        dec (float): Declination

    ReturnFs:
        x_component: the x components as an array where the first element corresponds 
        to the first set of ra and dec and continues in order
        y_component: the y components as an array where the first element corresponds 
        to the first set of ra and dec and continues in order
        z_component: the z components as an array where the first element corresponds 
        to the first set of ra and dec and continues in order
    """
    # flattening the variables and checking their sizes
    ra = flat_variable(ra)
    dec = flat_variable(dec)
    if len(ra) != len(dec):
        print("Error, calculate_unit_vector function on the file gauss_method.py")
        print("please, check the dimension of the values that you provided for the RA and DEC")
        print("ra values: ", ra)
        print("dec values: ", dec)
        raise ValueError
    # math
    x_component = np.cos(ra) * np.cos(dec) #x hat
    y_component = np.sin(ra) * np.cos(dec) #y hat
    z_component = np.sin(dec) #z hat
    return np.array([x_component,y_component,z_component])

def get_unit_matrix(ra,dec):
    """Provided RA and Dec values, this function separates them into sets of three in order. 
    If this is not possible, the program terminates. (The RA and Dec values may also already 
    be provided as sets of three.) For each set of three, it constructs a unit vector matrix 
    that represents the unit vectors of the given RA and Dec values.

    Args:
        ra (list): A list containing at least three float values representing the right ascension 
            of the object. The list must have the same length and shape as the argument `dec`.
        dec (list): A list containing at least three float values representing the declination 
            of the object. The list must have the same length and shape as the argument `ra`.

    Returns:
        np.ndarray: A matrix containing the unit vectors of the provided coordinates, each 
        set of three has corresponding a matrix in order.
    """
    # flattening the variables and checking their sizes
    ra = flat_variable(ra)
    dec = flat_variable(dec)
    if len(ra) != len(dec) or len(ra) % 3 != 0:
        print("Error, on function get_unit_matrix, file gauss_method.py")
        print("The provided RA and DEC values are not compatible")
        print("This are the provided RA values: ",ra)
        print("this are the provided DEC values: ",dec)
        print(len(ra),len(dec))
        raise ValueError
    
    x_cmp,y_cmp,z_cmp = unit_vector_from_ra_dec(ra,dec)
    unit_matrices = []
    
    # separating arrays every three elements
    x_cmp = x_cmp.reshape(-1,3)
    y_cmp = y_cmp.reshape(-1,3)
    z_cmp = z_cmp.reshape(-1,3)
    
    # checking the amount of create arrays are correct
    if len(x_cmp) != int(len(ra) / 3) and len(y_cmp) != int(len(ra) / 3) and len(z_cmp) != int(len(ra) / 3):
        print("Error, on function get_unit_matrix, file gauss_method.py")
        print("The provided RA and DEC values are not compatible")
        print("This are the provided RA values: ",ra)
        print("this are the provided DEC values: ",dec)
        print(x_cmp.shape)
        raise ValueError
        
    # Creating the matrices based on the separated arrays
    for i in range(int(len(ra) / 3)):
        current_matrix = np.array([x_cmp[i],y_cmp[i],z_cmp[i]])
        unit_matrices.append(current_matrix)
    return unit_matrices

def delta_t_values(t0,t1,t2):
    """Calculates the delta ts values of the gauss method

    Args:
        t0 (list / float): the first observation time
        t1 (list / float): the second observation time
        t2 (list / float): the third observation time

    Returns:
        tuple: a list of the tau values.
    """
    # flattening the provided input
    t0 = flat_variable(t0)
    t1 = flat_variable(t1)
    t2 = flat_variable(t2)
    if len(t0) != len(t1) or len(t0) != len(t2):
        print("Error, on function delta_t_values, file gauss_method.py")
        print("The provided values do not have the expected dimension")
        print(len(t0),len(t1),len(t2))
        raise ValueError
    # math
    d_t1 = t0 - t1 
    d_t3 = t2 - t1
    d_t = t2 - t0 
    return (d_t1,d_t3,d_t)

def calculate_a_and_b(d_t1,d_t3,d_t):
    """calculates a and b values, used on the coeffiecinets_for_polynomial function and on the gauss method.

    Args:
        d_t1 (list / float): delta one value
        d_t3 (list / float): delta three value
        d_t (list / float): tau or delta t value

    Returns:
        tuple: tuple with the a's and b's values.
    """
    # changing the given parameters into one dimensional (flat) numpy arrays
    d_t1 = flat_variable(d_t1)
    d_t3 = flat_variable(d_t3)
    d_t = flat_variable(d_t)
    # math
    a_1 = d_t3 / d_t
    b_1 = (d_t3 * (d_t**2 - d_t3**2)) / (6 * d_t)
    a_3 = - d_t1 / d_t
    b_3 = - (d_t1 * (d_t**2 - d_t1 **2)) / (6 * d_t)
    return (a_1,b_1,a_3,b_3)

def calculations_for_coefficients(a_1,b_1,a_3,b_3,position_R,rho,B_matrix):
    """calculate values used on the function coeffiecinets_for_polynomial.  
    
    Args:
        a_1 (list / float): a one value
        b_1 (list / float): b one value
        a_3 (list / float): a three value
        b_3 (list / float): b three value
        position_R (list / np.matrix): The position of the observer, 3x3
        rho (list / np.matrix): The unit vector of the matrix, 3x3
        B_matrix (list / np.matrix): matrix  used for the calculation, 3x3

    Returns:
        tuple: d_1 and d_2 values used to generate c values, and also additional 
        values used on automatic testing
    """
    a_1 = flat_variable(a_1)
    b_1 = flat_variable(b_1)
    a_3 = flat_variable(a_3)
    b_3 = flat_variable(b_3)
    # checking that the sizes of the parameters coincide 
    if len(a_1) != len(b_1) or len(a_1) != len(a_3) or len(a_1) != len(b_3):
        print("Error, on function calculations_for_coefficients, file gauss_method.py")
        print("The provided values do not have the expected dimension")
        print(len(a_1),len(b_1),len(a_3),len(b_3))
        raise ValueError
    # checking the dimension of the matrix parameters, and changing them if necessary    
    position_R = np.array(position_R)
    rho = np.array(rho)
    B_matrix = np.array(B_matrix)
    if position_R.ndim > 3:
        position_R = np.reshape(position_R,(len(a_1),3,3))
        rho = np.reshape(rho,(len(a_1),3,3))
        B_matrix = np.reshape(B_matrix,(len(a_1),3,3))
    if position_R.ndim == 2:
        position_R = np.reshape(position_R, (1,3,3))
        rho = np.reshape(rho, (1,3,3))
        B_matrix = np.reshape(B_matrix, (1,3,3))
    # separating necessary variables
    rh01 = np.array(rho[:,0,1])
    rh11 = np.array(rho[:,1,1])
    rh21 = np.array(rho[:,2,1])
    B_matrix = np.array(B_matrix,dtype = np.matrix)
    bmatrix10 = np.array(B_matrix[:,1,0])
    bmatrix11 = np.array(B_matrix[:,1,1])
    bmatrix12 = np.array(B_matrix[:,1,2])
    position_R = np.array(position_R,dtype = np.matrix)
    pos_r01 = np.array(position_R[:,0,1])
    pos_r11 = np.array(position_R[:,1,1])
    pos_r21 = np.array(position_R[:,2,1])
    temp_matrix1 = np.array([[rh01],[rh11],[rh21]],dtype = np.matrix).T # vector components of p_2
    temp_matrix2 = np.array([[pos_r01,pos_r11,pos_r21]],dtype = np.matrix).T # vector components of R_2
    # math
    d_1 = bmatrix10 * a_1 - bmatrix11 + bmatrix12 * a_3
    d_2 = bmatrix10 * b_1 + bmatrix12 * b_3
    part1_mag = pos_r01** 2 + pos_r11** 2 + pos_r21** 2
    part1_mag = part1_mag.astype(float)
    magnitude_R_2 = np.sqrt(part1_mag)
    pos_2_dot_p2 = temp_matrix1 @ temp_matrix2 
    pos_2_dot_p2 = pos_2_dot_p2.squeeze()
    return (d_1,d_2,magnitude_R_2,pos_2_dot_p2)

def coeffiecinets_for_polynomial(d_1,d_2,magnitude_R_2,pos_2_dot_p2):
    """Generates non-zero c values for the Gauss method c,c3,c6,c8 
    uses function calculations_for_coefficients. can be used for multiple sets 
    when parameters are provided as lit or stacks

    Args:
        a_1 (list / float): a one value
        b_1 (list / float): b one value
        a_3 (list / float): a three value
        b_3 (list float): b three value
        position_R (list / np.matrix): The position of the observer, 3x3 matrix
        rho (list / np.matrix): The unit vector of the matrix, 3x3 matrix
        B_matrix (list / np.matrix): matrix  used for the calculation 3x3 matrix
        
    Returns:
        tuple: with generated c values in the order c,c3,c6,c8
    """
    # flattening the provided parameters that could be given as a list
    d_1 = flat_variable(d_1)
    d_2 =  flat_variable(d_2)
    magnitude_R_2 = flat_variable(magnitude_R_2)
    pos_2_dot_p2 = flat_variable(pos_2_dot_p2)
    
    # checking that the size of the provided values coincide
    if len(d_1) != len(d_2) or len(d_1) != len(magnitude_R_2) or len(d_1) != len(pos_2_dot_p2):
        print("Error, on function coeffiecinets_for_polynomial, file gauss_method.py")
        print("The provided values do not have the expected dimension")
        print(len(d_1),len(d_2),len(magnitude_R_2),len(pos_2_dot_p2))
        raise ValueError
    
    # math
    c = d_2 ** 2 * mue**2
    c_3 =mue * 2 * (d_1 * d_2 + d_2 * pos_2_dot_p2)
    c_6 = d_1**2 +magnitude_R_2**2 + 2 * d_1 *pos_2_dot_p2
    c_8 = np.full((len(d_1)), -1)    # it can be -1 and, or multiply all the other coefficients by -1 instead.
    c = c.squeeze()
    c_3 = c_3.squeeze()
    c_6 = c_6.squeeze()
    c_8 = c_8.squeeze()
    return (c,c_3,c_6,c_8)

def rho_magnitude_and_r_values(root,B_matrix,a_1,b_1,a_3,b_3,R,p):
    """This function calculates the rho and r values of the object. 
    it can handle multiple objects in a single call, when parameters 
    are given as list and stacks
    
    Args:
        root (list / float): calculated using the current root in question (r_2)
        B_matrix (list / np.array): matrix  used for the calculation
        a_1 (list / float): a one value
        b_1 (list / float): b one value
        a_3 (list / float): a three value
        b_3 (list / float): b three value
        R (list / np.matrix): The position of the observer 
        P (list / np.matrix): The unit vector of the matrix
        
    Returns:
        list: list with 3 of lenght 
            * The firt element is the root, 
            * The second element contains all the calculated rho values (magnitudes),
            as an list of leght three.
            * The third element contains calculated (r) position values of the object 
            as a value is a 3x3 matrix per each set 
            (if multiple sets are present the matrixes are stacked in order)
    """
    # defining giving terms as numpy arrays and flattening them
    root = flat_variable(root)
    a_1 = flat_variable(a_1)
    b_1 = flat_variable(b_1)
    a_3 = flat_variable(a_3)
    b_3 = flat_variable(b_3)
    B_matrix = np.array(B_matrix)
    p = np.array(p)
    R = np.array(R)
    # cheking that the lenght of the provided values is the same
    if len(root) != len(a_1) or len(a_1) != len(b_1) or len(a_1) != len(a_3) or len(a_1) != len(b_3):
        print("Error, on function calculations_for_coefficients, file gauss_method.py")
        print("The provided values do not have the expected dimension")
        print(len(root),len(a_1),len(b_1),len(a_3),len(b_3))
        raise ValueError  
    
    # cheking and changing diminsion when neccesary
    if p.ndim > 3:
        p = np.reshape(p,(len(a_1),3,3))
        R = np.reshape(R,(len(a_1),3,3))
        B_matrix = np.reshape(B_matrix,(len(a_1),3,3))
    if p.ndim == 2:
        p = np.reshape(p, (1,3,3))
        R = np.reshape(R, (1,3,3))
        B_matrix = np.reshape(B_matrix, (1,3,3))
        
    # defining terms that will be used
    sigma = np.array(mue/(root**3),dtype=float).squeeze()
    b_matrix00 = B_matrix[:,0,0]
    b_matrix10 = B_matrix[:,1,0]
    b_matrix20 = B_matrix[:,2,0]
    b_matrix01 = B_matrix[:,0,1]
    b_matrix11 = B_matrix[:,1,1]
    b_matrix21 = B_matrix[:,2,1]
    b_matrix02 = B_matrix[:,0,2]
    b_matrix12 = B_matrix[:,1,2]
    b_matrix22 = B_matrix[:,2,2]
    r0 = R[:,:,0]
    r1 = R[:,:,1]
    r2 = R[:,:,2]
    p0 = p[:,:,0]
    p1 = p[:,:,1]
    p2 = p[:,:,2]
    r = np.zeros((len(a_1),3,3))
    
    # math
    mag_p1 = (b_matrix00 * a_1 - b_matrix01 + b_matrix02 * a_3 + (b_matrix00 * b_1 + b_matrix02 * b_3) * sigma)/((-1)*(a_1+b_1*sigma))
    mag_p2 = (b_matrix10 * a_1 - b_matrix11 + b_matrix12 * a_3 + (b_matrix10 * b_1 + b_matrix12 * b_3) * sigma)
    mag_p3 = (b_matrix20 * a_1 - b_matrix21 + b_matrix22 * a_3 + (b_matrix20 * b_1 + b_matrix22 * b_3) * sigma)/ ((-1)* (a_3+b_3*sigma))
    r[:,:,0] = np.array(r0) + np.array(np.matrix(mag_p1).T)*np.array(p0)
    r[:,:,1] = np.array(r1) + np.array(np.matrix(mag_p2).T)*np.array(p1)
    r[:,:,2] = np.array(r2) + np.array(np.matrix(mag_p3).T)*np.array(p2)
    magnitudes = np.array([mag_p1.squeeze(),mag_p2.squeeze(),mag_p3.squeeze()])
    return (root, magnitudes,r)

def Check_root(rho2,spurious_distance= 0.001):
    """Checks for SPurious roots with Spurious distance being the 
    earth Hill's sphere

    Args:
        rho2 (list / float): The rho2 generated by the root that is being tested. 
        The parameter can be proviede as a float the object
        spurious_distance (float, optional): the Eartch Hill's sphere in AU. 
        Defaults to 0.001.

    Returns:
        list: a list of booleans values where each rerpesnt each of the root 
        tested, the roots that correspond to the True value spurious root
        """
    rho2 = flat_variable(rho2)
    bool_value = rho2 < spurious_distance # staments that makes the check 
    return bool_value

def provide_position_matrix():
    """Functions where user can provided the position matric of the obser 
    by typing it out, multiple matrices can be written as a stakc. 

    Returns:
        list: unit matrices stored within the functio
    """
    postion_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
    postion_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
    postion_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
    postion_matrix = np.array([[postion_vector1[0],postion_vector2[0],postion_vector3[0]],
                            [postion_vector1[1],postion_vector2[1],postion_vector3[1]],
                            [postion_vector1[2],postion_vector2[2],postion_vector3[2]]])
    return postion_matrix

def provide_unit_matrix():
    """Functions where user can provided the unit matrix of the observation by typing it out. 

    Returns:
        list: unit matrices stored within the functio
    """
    unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
    unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
    unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
    unit_matrix = np.array([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],
                            [unit_vector1[1],unit_vector2[1],unit_vector3[1]],
                            [unit_vector1[2],unit_vector2[2],unit_vector3[2]]])
    return unit_matrix

def calculated_v2(position_matrix,t1,t3):
    """This is a function to valculate the velocity vector v2 of an object.
    it can perform the calculations for multiple object if parameters 
    are given as list, it uses gibbs method to calculate the v2 values

    Args:
        position_matrix (list): a three by three matrix with the position of the object. 
        it can be given as a stack of 3x3 matrices to perform calculations in multiple sets
        t1 (list / float): The first time of observation. If calculations are perform in
        multiple sets the parameter is a list of floats
        t3 (list / float): The third time of observation. If calculations are perform in 
        multiple sets the parameter is a list of floats

    Returns:
        list: a list of lenght of three that represent the v2 vector for the object(s)
    """
    # Converting provided values to numpy arryas
    t1 = np.array(t1)
    t1 = np.atleast_1d([t1]) # Result keep [] brackets and changed the shape of it
    t1 = t1.flatten()
    t3 = np.array(t3)
    t3 = np.atleast_1d([t3]) # Result keep [] brackets and changed the shape of it
    t3 = t3.flatten()
    position_matrix = np.array(position_matrix)
    
    # cheking that the lenght of the provided values are the same
    if len(t1) != len(t3):
        print("Error, on function calculate_v2, file gauss_method.py")
        print("The provided values do not have the expected dimension")
        print(len(t1),len(t3))
        raise ValueError
        
    # cheking and changing diminsion if neccesary to have a stack of matrices
    if position_matrix.ndim > 3:
        position_matrix = np.reshape(position_matrix,(len(t1),3,3))
        
    if position_matrix.ndim == 2:
        position_matrix = np.reshape(position_matrix, (1,3,3))
    
    # math
    t_1 = np.sqrt(mue) * t1
    t_3 = np.sqrt(mue) * t3
    t13 = t_3 - t_1
    norm = []
    for matrix in position_matrix:
        matrix = np.array(matrix,dtype=float)
        norm.append( np.array(1/(np.array([np.linalg.norm(matrix,axis=0)]))**3).squeeze())
    norm = np.array(norm)
    r1m3 = norm[:,0]
    r2m3 = norm[:,1]
    r3m3 = norm[:,2]
    d1= t_3 * (r1m3 / 12 - 1/ (t_1 * t13)) 
    d2= (t_1+t_3) * (r2m3 / 12 - 1/ (t_1 * t_3)) 
    d3= -t_1 * (r3m3 / 12 + 1/ (t_3 * t13)) 
    v2 =[]
    for i in range(3) :
        v2.append(np.sqrt(mue)*(-d1*position_matrix[:,i,0]+d2*position_matrix[:,i,1]+d3*position_matrix[:,i,2]))
    v2 = np.array(v2).T
    v2 = np.array(v2).squeeze()
    return v2

def get_roots_values(roots,B_matrix,delta_t1,delta_t3,a1,b1,a3,b3,observer_pos,unit_vector_rho):
    """This function calculates various values based on input parameters 
    and returns them, it utuilizes the functions, rho_magnitude_and_r_values to calculate the
    rho magnitudes and r values of the root, function calculated_v2 to calculated the 
    velocity at time two of the root, and the function check_root, too see if the 
    roots provided are spurious.

    Args:
        roots (list): a list of the roots
        B_matrix (np.array): the array containg the B matrix or /matrices 
        if we are doing multiple sets
        delta_t1 (int/list): delta_t1 values
        delta_t3 (int/list): delta_t3 values
        a1 (int/list): a1 value
        b1 (int/list): b1 value
        a3 (int/list): a3 value
        b3 (int/list): b3 value
        observer_pos (np.array): the array conataining the observer position
        unit_vector_rho (np.array): the array conataining the unit vector of rho

    Returns:
        tuple: containing the following values:
            - test_root: a spurious root check returning a true or false array, 
            where true means that the root is spurious
            - r_values: containg the values calculated form function rho_magnitude_and_r_values
            - v2: the velocity of the asteroid at the second observation time, calculated 
            using gibbs method and function calculated_v2
            - r_2s: the r(position of the asteroid from the center of the sun) value at 
            the second observation time.
    """
    new_B_matrix = np.repeat(B_matrix[:, None], 3, axis=1)
    new_a1 = np.repeat(a1[:, None], 3, axis=1)
    new_b1 = np.repeat(b1[:, None], 3, axis=1)
    new_a3 = np.repeat(a3[:, None], 3, axis=1)
    new_b3 = np.repeat(b3[:, None], 3, axis=1)
    new_R = np.repeat(observer_pos[:, None], 3, axis=1)
    new_p = np.repeat(unit_vector_rho[:, None], 3, axis=1)
    current_roots,magnitudes,r_values = rho_magnitude_and_r_values(roots,new_B_matrix,new_a1,new_b1,new_a3,new_b3,new_R,new_p)
    v2_delta_t1 = np.repeat(delta_t1[:, None], 3,axis=1)
    v2_delta_t3 = np.repeat(delta_t3[:, None], 3,axis=1) 
    v2 = calculated_v2(r_values,v2_delta_t1,v2_delta_t3)
    r_2s = r_values[:,:,1]
    rho_2s = magnitudes[1]
    test_root = Check_root(rho_2s) # chedking for spurious roots
    return test_root,r_values,v2,r_2s

def gauss_method(observation_times,unit_vector_rho,observer_position):
    """Performs the Gauss method on either a single set or multiple sets of data.  
    This function relies on the 'math_functions' file, and on local functions  
    'delta_t_values', 'calculate_a_and_b', 'coeffiecinets_for_polynomial', and 'rho_magnitude_and_r_values'.  

    Args:  
        observation_times (list): A list of modified Julian dates. The total number of dates  
            must be divisible by three, since they are grouped into sets of three in the  
            order provided.  
        unit_vector (np.ndarray): The unit vector matrix of the asteroid.  
        observer_position (np.ndarray): The position matrix of the observatory.  
        
    Returns:
        list: A list containing the results of the Gauss method with the following elements:
            
            - **roots** (list): A list of the calculated roots in order 
            (e.g., [root1, root2, root3, ...]).
            
            - **test_root** (list of bool): A list indicating whether each root is spurious. 
            For example: [True, False, False], where `True` marks a spurious root.
            
            - **r_values** (list of arrays): If exactly three roots are found, this contains 
            three 3x3 matrices representing the position vectors of the asteroid.
        
            - **observation_times** (list of floats): containst the tisl of observation times.
            
            - **v2** (list of arrays): A list of velocity vectors (length-3 arrays). 
            If three roots are found, the shape will be (3, 3). 
            
            - **r2_s** (list of arrays): Orbital elements for all roots, structured 
            as separate arrays for each element type in the order 
            [a, e, i, long_node, arg_peri, mean_anomaly].  
            Each array contains the values for all roots in order.  
            Example: [[a1, a2, a3], [e1, e2, e3], [i1, i2, i3], ...].
    """
    
    # falteening the provided observatory times, and cheking the the lenght provided inputs match
    size_tobs = flat_variable(observation_times)
    if len(size_tobs) % 3 != 0 or len(observer_position) != len(unit_vector_rho) or len(observer_position) != len(observation_times):
        print("Error, on function gauss_method, file gauss_method.py")
        print("The provided values do not have the expected dimension")
        print(len(observation_times),len(observer_position),len(unit_vector_rho))
        raise ValueError
    
    # size determines the amount of sets we are working with
    size = int(len(size_tobs) / 3)
    
    # making R, and P numpy arrays and changing dimensions if neccerary
    unit_vector_rho = np.array(unit_vector_rho)
    observer_position = np.array(observer_position)
    if observer_position.ndim > 3:
        observer_position = np.reshape(observer_position,(size,3,3))
        unit_vector_rho = np.reshape(unit_vector_rho,(size,3,3))
    if observer_position.ndim == 2:
        observer_position = np.reshape(observer_position, (1,3,3))
        unit_vector_rho = np.reshape(unit_vector_rho, (1,3,3))
    observation_times = np.array(observation_times)
    
    # dividing observation times into set of three
    observation_times = np.reshape(observation_times, (size,3))
    
    # for each set of three separating the values into t0 t1 and t2.
    t0 = observation_times[:,0]
    t1 = observation_times[:,1]
    t2 = observation_times[:,2]
    

    # math and function calls
    inverse_unit_rho = mf.take_inverse_matrix(unit_vector_rho)
    B_matrix = inverse_unit_rho@observer_position
    delta_t1,delta_t3,delta_t = delta_t_values(t0,t1,t2)
    a1,b1,a3,b3 = calculate_a_and_b(delta_t1,delta_t3,delta_t)
    d_1,d_2,magnitude_R_2,pos_2_dot_p2 = calculations_for_coefficients(a1,b1,a3,b3,observer_position,unit_vector_rho,B_matrix)
    c,c3,c6,c8 = coeffiecinets_for_polynomial(d_1,d_2,magnitude_R_2,pos_2_dot_p2)
    roots = mf.eight_equation_four_coeff(c,c3,c6,c8)
    
    # getting the important values for each of the three positive real roots generated. 
    test_root,r_values,v2,r_2s = get_roots_values(roots,B_matrix,delta_t1,delta_t3,a1,b1,a3,b3,observer_position,unit_vector_rho)
    observation_times = observation_times.squeeze()
    return roots,test_root,r_values,observation_times,v2,r_2s

def run_gauss_method_hand_values(date_time,obscode=500):
    """Given the position values, unit vector, observation times, and observatory code,  
    this function calculates the Gauss method solution and the orbital elements.  
    It makes use of the helper functions `provide_unit_matrix` and `provide_position_matrix`,  
    which are expected to have already contain correct matrices.
    
    Args:
        date_time (list): A list of observation times. At least three observation times 
            must be provided. For Gauss' method, the times are grouped into sets of three 
            in the order they are given.
        obscode (str): Observatory code (only one observatory code can be provided), defaults
            to 500. and is not used for the calulations.

    Returns:
        list: A list containing the results of the Gauss method with the following elements:
            
            - **roots** (list): A list of the calculated roots in order 
            (e.g., [root1, root2, root3, ...]).
            
            - **test_root** (list of bool): A list indicating whether each root is spurious. 
            For example: [True, False, False], where `True` marks a spurious root.
            
            - **r_values** (list of arrays): If exactly three roots are found, this contains 
            three 3x3 matrices representing the position vectors of the asteroid.
            
            - **v2** (list of arrays): A list of velocity vectors (length-3 arrays). 
            If three roots are found, the shape will be (3, 3).
            
            - **orbit_elements** (list of arrays): Orbital elements for all roots, structured 
            as separate arrays for each element type in the order 
            [a, e, i, long_node, arg_peri, mean_anomaly].  
            Each array contains the values for all roots in order.  
            Example: [[a1, a2, a3], [e1, e2, e3], [i1, i2, i3], ...].
    """
    postion_matrix = provide_position_matrix() 
    unit_matrix = provide_unit_matrix()
    result = gauss_method([date_time,date_time],[unit_matrix,unit_matrix],[postion_matrix,postion_matrix])
    return result

def run_gauss_method_wis(date_time,ra_dec_values,obscode):
    """ Applies Gauss' method using the functions 'calculate_position_obs', 'get_unit_matrix', 
    'gauss_method', and the library 'wis' to estimate orbital elements from three or more 
    astrometric observations.

    Args:
        date_time (list): A list of observation times. At least three observation times 
            must be provided. For Gauss' method, the times are grouped into sets of three 
            in the order they are given.
        ra_dec_values (list): A list of sublists containing right ascension and declination 
            values in the format [RA, DEC]. Each set of three times of observation must have a 
            corresponding [RA, DEC] entry. 
        obscode (str): Observatory code (only one observatory code can be provided).

    Returns:
        list: A list containing the results of the Gauss method with the following elements:
            
            - **roots** (list): A list of the calculated roots in order 
            (e.g., [root1, root2, root3, ...]).
            
            - **test_root** (list of bool): A list indicating whether each root is spurious. 
            For example: [True, False, False], where `True` marks a spurious root.
            
            - **r_values** (list of arrays): If exactly three roots are found, this contains 
            three 3x3 matrices representing the position vectors of the asteroid.
            
            - **v2** (list of arrays): A list of velocity vectors (length-3 arrays). 
            If three roots are found, the shape will be (3, 3).
            
            - **orbit_elements** (list of arrays): Orbital elements for all roots, structured 
            as separate arrays for each element type in the order 
            [a, e, i, long_node, arg_peri, mean_anomaly].  
            Each array contains the values for all roots in order.  
            Example: [[a1, a2, a3], [e1, e2, e3], [i1, i2, i3], ...].
    """
    postion_matrix = position_obs_wis(date_time.copy(),obscode) 
    postion_matrix = postion_matrix.T 
    unit_matrix = get_unit_matrix(ra_dec_values) # [[alpha,delta][]...]
    result = gauss_method(date_time,unit_matrix,postion_matrix)  
    return result
