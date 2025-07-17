# Description: Contains math functions that implement numpy such as dot products, matrix multiplication, and matrix inverses. 

import numpy as np


def take_inverse_matrix(matrix):
    try:
        inverse_matrix = np.linalg.inv(matrix)
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
    try:
        mult = np.dot(matrix1,matrix2)
    except:
        print(f"Error: The matrix multiplication was not possible please double check the matrices dimension and values {matrix1}{matrix2}")
        print("function: matrix_dot_prod")
        exit()
    return mult

def eight_equation_four_coeff(c,c3,c6,c8):
    c1 = 0
    c2 = 0
    c4 =0
    c5 = 0
    c7 = 0
    coefficients = [float(c8),c7,float(c6),c5,c4,float(c3),c2,c1,float(c)] # some numbers are passed in as matrices elements thus making them floats.
    try:
        roots = np.roots(coefficients)
    except:
        print(f'Unable to produce a root: {coefficients}')
        print("function test_four_coeff")
        exit()
    real_roots = []
    roots_count = []
    positive_roots = []
    for root in roots:
        if np.isreal(root):
            real_roots.append(float(root))
            roots_count.append(True)
            if float(root) >= 0.0:
                positive_roots.append(float(root))
        else:
            roots_count.append(False)
    return positive_roots

