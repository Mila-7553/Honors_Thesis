'''This file contains functions that make use of math utility functions built on NumPy, 
including dot products, matrix multiplications, and matrix inverses.
'''
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=np.exceptions.ComplexWarning)

def take_inverse_matrix(matrix):
    """This function uses numpy.linalg.inv and numpy.linalg.pinv to take the inverse of a provided matrix

    Args:
        matrix (np.ndarray): a matrix to be inversed

    Returns:
        np.ndarray: the inverse matrix
    """
    try:
        inverse_matrix = np.linalg.inv(matrix)
    except:
        try:
            inverse_matrix = np.linalg.pinv(matrix) # pseudo inverse function by numpy
        except (np.linalg.LinAlgError,ValueError) as e:
            if isinstance(e, np.linalg.LinAlgError):
                print(f"Error: the matrix provided is not invertible {matrix}")
                print("function: take_inverse")
            else:
                print(f'Error: the provided data {matrix} is not a matrix please provide different data.')
                print("function: take_inverse")
            exit() # can change to return a flag
    return inverse_matrix

def matrix_dot_prod(matrix1,matrix2):
    """Function to take the matrix dot product of the provided matrix1 dot provided matrix2

    Args:
        matrix1 (np.matrix): matrix to be dotted
        matrix2 (np.matrix): matrix to be dotted

    Returns:
        np.ndarray: result of the dot product
    """
    # only used by fortran replica
    try:
        mult = np.dot(matrix1,matrix2)
    except:
        print(f"Error: The matrix multiplication was not possible please double check the matrices dimension and values {matrix1}{matrix2}")
        print("function: matrix_dot_prod")
        exit()
    return mult

def eight_equation_four_coeff(c,c3,c6,c8):
    """Function to get the roots for a polynomial of 8 degree provided the coefficient 0, 3, 6 and 8

    Args:
        c (np.array or float): the provided coefficient c
        c3 (np.array or float): the provided coefficient c3
        c6 (np.array or float): the provided coefficient c6
        c8 (np.array or float): the provided coefficient c8

    Returns:
        np.array: returns the positive and real roots generated from the provided coefficients
    """
    c = np.atleast_1d(c)
    c = c.flatten()
    c3 = np.atleast_1d(c3)
    c3 = c3.flatten()
    c6 = np.atleast_1d(c6)
    c6 = c6.flatten()
    c8 = np.atleast_1d(c8)
    c8 = c8.flatten()
    if len(c) != len(c3) or len(c) != len(c6) or len(c) != len(c8):
        print("Error, on function eight_equation_four_coeff, file gauss_method.py")
        print("The provided values do not have the expected dimension")
        print(len(c),len(c3),len(c6),len(c8))
        raise ValueError
    size = len(c)
    if type(c) == list or type(c) == np.ndarray and len(c6) == len(c8):
        size = len(c6)
    else:
        size = 1
    # print(c)
    c = np.array(c,dtype=np.float64).squeeze()
    c3 = np.array(c3,dtype=np.float64).squeeze()
    c6 = np.array(c6,dtype=np.float64).squeeze()
    c8 = np.array(c8,dtype=np.float64).squeeze()
    c1 = np.full(size, 0).squeeze()
    c2 = np.full(size, 0).squeeze()
    c4 = np.full(size, 0).squeeze()
    c5 = np.full(size, 0).squeeze()
    c7 = np.full(size, 0).squeeze()
    coefficients = np.array([c8,c7,c6,c5,c4,c3,c2,c1,c]) # some numbers are passed in as matrices elements thus making them floats.
    coefficients = coefficients.T
    if size == 1:
        coefficients = [coefficients]
    try:
        roots_all = np.array([np.roots(cs) for cs in coefficients], dtype=complex)
    except:
        print(f'Unable to produce a root: {coefficients}')
        print("function test_four_coeff")
        exit()
    positive_roots = [] # still need to test this values are able to run without errors. 
    # positive_roots = [-100,-100,-100] # still need to test this values are able to run without errors. 
    # positive_roots = [roots[np.isreal(roots) & (roots.real >= 0)].real
    # for roots in roots_all] # makes sure they stay in their apporpiate sublist
    # positive_roots = np.array(positive_roots,dtype=np.longdouble).squeeze()

    for i in range(len(roots_all)):
        curren_pos = roots_all[i][np.isreal(roots_all[i]) & (roots_all[i] > 0)]
        while len(curren_pos) < 3:
            curren_pos = np.append(curren_pos,0.0001)
        curren_pos = curren_pos.astype(np.float64)
        positive_roots.append([curren_pos])
    positive_roots = np.array(positive_roots)
    positive_roots = positive_roots.squeeze()
    return positive_roots
