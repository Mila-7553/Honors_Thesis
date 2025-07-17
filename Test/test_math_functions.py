import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../Honor_thesis_Code")))
import math_fucntions as mf




unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
test_matrix1 = np.matrix([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],[unit_vector1[1],unit_vector2[1],unit_vector3[1]],[unit_vector1[2],unit_vector2[2],unit_vector3[2]]])
# should have a t least three

def test_take_inverse_matrix():
    results = mf.take_inverse_matrix(test_matrix1)
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
    
# test_take_inverse_matrix()

#  Fortran results:
#   Edit: inverse unit Vector 1:   -181.39997132286484        346.82086618012170       -167.22601164748204     
#   Edit: inverse unit Vector 2:   -49.141453222161964        99.181666752383251       -49.583095210648771     
#   Edit: invers unit Vector 3:    732.74446483456154       -1365.9487590320009        637.88133606940517 

# Python results:
# -181.39997132286484   346.8208661801217   -167.22601164748204
# -49.141453222161964   99.18166675238325   -49.58309521064877
# 732.7444648345615   -1365.948759032001   637.8813360694052
# not as precise example at [2,2]



matrix1_test1 = [[8,4,5],[6,3,7],[5,7,3]]
matrix2_test1 = [[5],[3],[1]]
matrix1_test2 = [[1,0,0],[0,1,0],[0,0,1]]
matrix2_test2 = 5

postion_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
postion_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
postion_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
matrix1_test3 = [[postion_vector1[0],postion_vector2[0],postion_vector3[0]],[postion_vector1[1],postion_vector2[1],postion_vector3[1]],[postion_vector1[2],postion_vector2[2],postion_vector3[2]]]
matrix2_test3 = [[0.0032114744736032957],[0],[0.0030602457022409226]]

matrix1_test4 = matrix1_test3.copy()
matrix2_test4 = [[0.5361692088746649],[-1],[0.46383079112533504]]

def test_matrix_dot_prod():
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
    
    result3 = mf.matrix_dot_prod(matrix1_test3,matrix2_test3)
    assert len(result3) == 3
    assert float(result3[0]) == -5.8515872213727631E-003
    assert float(result3[1]) == -1.9995299459034648E-003
    assert float(result3[2]) == -8.6677509347662051E-004
    
    result4 = mf.matrix_dot_prod(matrix1_test4,matrix2_test4)
    assert len(result4) == 3
    assert float(result4[0]) == 5.8502493672572542E-003
    assert float(result4[1]) == 2.0256251869944231E-003
    assert float(result4[2]) == 8.7695014435107410E-004
# test_matrix_dot_prod()
# fortran results:
#  Edit rb:   -5.8515872213727631E-003  -1.9995299459034648E-003  -8.6677509347662051E-004

# python results:
# -0.005851587221372763
# -0.001999529945903465    !
# -0.0008667750934766205

# fortran results
#  Edit ra:    5.8502493672572542E-003   2.0256251869944231E-003   8.7695014435107410E-004

# python results
# 0.005850249367257254
# 0.002025625186994423
# 0.0008769501443510741
# Not as precise (does not include last digit)



test1_c = -1.0895161396889805
test1_c3 =2.8036979214270064
test1_c6 =-2.7106217754654516
test1_c8 =1
def test_eight_equation_four_coeff():
    result1 = mf.eight_equation_four_coeff(test1_c,test1_c3,test1_c6,test1_c8)
    assert round(float(result1[0]),13) == round(1.2457493693119874,13)
    assert round(float(result1[1]),13) == round(1.0544796939362728,13)
    assert round(float(result1[2]),13) == round(0.98366052181486430,13)

# fortran roots:
# 0.98366052181486430   (spurious)
# 1.0544796939362728
# 1.2457493693119874

# python solutions
# 0.9836605218148519
# 1.054479693936294
# 1.245749369311985
