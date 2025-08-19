'''Contains math functions that implement numpy such as dot products, matrix multiplication, and matrix inverses. 
'''

import numpy as np
import warnings
warnings.filterwarnings("ignore", category=np.exceptions.ComplexWarning)


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
    # only used by fortran replica
    try:
        mult = np.dot(matrix1,matrix2)
    except:
        print(f"Error: The matrix multiplication was not possible please double check the matrices dimension and values {matrix1}{matrix2}")
        print("function: matrix_dot_prod")
        exit()
    return mult

def eight_equation_four_coeff(c,c3,c6,c8):
    if type(c) == list or type(c) == np.ndarray and len(c6) == len(c8):
        size = len(c6)
    else:
        size = 1
    c = np.array(c,dtype=np.float64)
    c3 = np.array(c3,dtype=np.float64)
    c6 = np.array(c6,dtype=np.float64)
    c8 = np.array(c8,dtype=np.float64)
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
    positive_roots = [-100,-100,-100] # still need to test this values are able to run without errors. 
    positive_roots = [roots[np.isreal(roots) & (roots.real >= 0)].real
    for roots in roots_all] # makes sure they stay in their apporpiate sublist
    return positive_roots

onec = 1.0895161396889808
onec3 = -2.8036979214269007
onec6 = 2.7106217754653157
onec8 = -1
# eight_equation_four_coeff(onec,onec3,onec6,onec8)

twoc = 10.940364370765677
twoc3 = -16.968656069203668
twoc6 = 6.722062703325108
twoc8 = -1

# eight_equation_four_coeff([onec,twoc],[onec3,twoc3],[onec6,twoc6],[onec8,twoc8])