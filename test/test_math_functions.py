import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../Honor_thesis_Code")))
from gauss import math_fucntions as mf
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
unit_matrix = np.matrix([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],
                        [unit_vector1[1],unit_vector2[1],unit_vector3[1]],
                        [unit_vector1[2],unit_vector2[2],unit_vector3[2]]])
matrix1_test1 = [[8,4,5],[6,3,7],[5,7,3]]
matrix2_test1 = [[5],[3],[1]]
matrix1_test2 = [[1,0,0],[0,1,0],[0,0,1]]
matrix2_test2 = 5
position_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
position_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
position_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
position_matrix = [[position_vector1[0],position_vector2[0],position_vector3[0]],
                [position_vector1[1],position_vector2[1],position_vector3[1]],
                [position_vector1[2],position_vector2[2],position_vector3[2]]]
matrix2_test3 = [[0.0032114744736032957],[0],[0.0030602457022409226]]
matrix2_test4 = [[0.5361692088746649],[-1],[0.46383079112533504]]
test1_c = -1.0895161396889805
test1_c3 =2.8036979214270064
test1_c6 =-2.7106217754654516
test1_c8 =1

def test_take_inverse_matrix():
    """A test for the function take_inverse_matrix on the file math_functions.py, that provided
    two matrices. test are performed for the return value of the function. 
    
    Provided values
        unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
        unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
        unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
        unit_matrix = np.matrix([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],
                                [unit_vector1[1],unit_vector2[1],unit_vector3[1]],
                                [unit_vector1[2],unit_vector2[2],unit_vector3[2]]])
            
    Results:
        The code result are ->
        inverse unit Vector 1 = -181.39997132286484   346.8208661801217   -167.22601164748204
        inverse unit Vector 2 = -49.141453222161964   99.18166675238325   -49.58309521064877
        inverse unit Vector 3 = 732.7444648345615   -1365.948759032001   637.8813360694052

        The fortran results are ->
        inverse unit Vector 1:   -181.39997132286484        346.82086618012170       -167.22601164748204     
        inverse unit Vector 2:   -49.141453222161964        99.181666752383251       -49.583095210648771     
        inverse unit Vector 3:    732.74446483456154       -1365.9487590320009        637.88133606940517 
    """
    results = mf.take_inverse_matrix(unit_matrix)
    assert type(results) == np.matrix
    assert len(results) == 3
    assert len(results[:]) == 3
    assert float(results[0,0]) == -181.39997132286484
    assert float(results[1,0]) == 346.82086618012170
    assert float(results[2,0]) == -167.22601164748204
    assert float(results[0,1]) == -49.141453222161964
    assert float(results[1,1]) == 99.181666752383251
    assert float(results[2,1]) == -49.583095210648771
    assert float(results[0,2]) == 732.74446483456154
    assert float(results[1,2]) == -1365.9487590320009
    assert float(results[2,2]) == 637.88133606940517

def test_matrix_dot_prod():
    """A test for the function matrix_dot_prod on the file math_functions.py, that provided
    two matrices. test are performed for the return value of the function. 
    
    Provided values
        matrix1_test1 = [[8,4,5],[6,3,7],[5,7,3]]
        matrix2_test1 = [[5],[3],[1]]
        matrix1_test2 = [[1,0,0],[0,1,0],[0,0,1]]
        matrix2_test2 = 5
        postion_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
        postion_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
        postion_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
        position_matrix = [[postion_vector1[0],postion_vector2[0],postion_vector3[0]],
        [postion_vector1[1],postion_vector2[1],postion_vector3[1]],[postion_vector1[2],
        postion_vector2[2],postion_vector3[2]]]
        matrix2_test3 = [[0.0032114744736032957],[0],[0.0030602457022409226]]
        matrix2_test4 = [[0.5361692088746649],[-1],[0.46383079112533504]]
    
    Results:
        The code result are ->
        0.005850249367257254
        0.002025625186994423
        0.0008769501443510741
                
        The fortran results are ->
        5.8502493672572542E-003   
        2.0256251869944231E-003   
        8.7695014435107410E-004
    """
    result1 = mf.matrix_dot_prod(matrix1_test1,matrix2_test1)
    assert len(result1) == 3
    assert float(result1[0,0]) == 57
    assert float(result1[1,0]) == 46
    assert float(result1[2,0]) == 49
    
    result2 = mf.matrix_dot_prod(matrix1_test2,matrix2_test2)
    assert len(result2) == 3
    assert len(result2[:]) == 3
    assert float(result2[0,0]) == 5
    assert float(result2[1,1]) == 5
    assert float(result2[2,2]) == 5
    
    result3 = mf.matrix_dot_prod(position_matrix,matrix2_test3)
    assert len(result3) == 3
    assert float(result3[0]) == -5.8515872213727631E-003
    assert float(result3[1]) == -1.9995299459034648E-003
    assert float(result3[2]) == -8.6677509347662051E-004
    
    result4 = mf.matrix_dot_prod(position_matrix,matrix2_test4)
    assert len(result4) == 3
    assert float(result4[0]) == 5.8502493672572542E-003
    assert float(result4[1]) == 2.0256251869944231E-003
    assert float(result4[2]) == 8.7695014435107410E-004

def test_eight_equation_four_coeff():
    """A test for the function eight_equation_four_coeff on the file math_functions.py, that provided test 
    coefficients c, c3, c6 and c8. test are performed for the return value of the function. 
    
    Provided values
        test1_c = -1.0895161396889805
        test1_c3 =2.8036979214270064
        test1_c6 =-2.7106217754654516
        test1_c8 =1
    
    Results:
        The code result are ->
        0.9836605218148519
        1.054479693936294
        1.245749369311985
        
        The fortran results are ->
        0.98366052181486430   
        1.0544796939362728
        1.2457493693119874
    """
    result1 = mf.eight_equation_four_coeff(test1_c,test1_c3,test1_c6,test1_c8)
    assert round(float(result1[0]),13) == round(1.2457493693119874,13)
    assert round(float(result1[1]),13) == round(1.0544796939362728,13)
    assert round(float(result1[2]),13) == round(0.98366052181486430,13)
