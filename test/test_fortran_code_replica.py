'''
Test script for fortran_code_replica.py.

This file runs tests using only test sample 1, and the results from the python code. 
The calculations are performed entirely in Python (no Fortran execution).
It also includes a comparison between the Python-generated results 
and the corresponding Fortran results for validation.
The fortran results were obtained by running the orbfit fortran code with the test sample 1 data.
'''
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../Honor_thesis_Code")))
from gauss import fortran_code_replica as fcr
import numpy as np

# [Add a description of what these unit_vector* quantities are and how they were obtained ... ]
unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
unit_vector_matrix1 = np.matrix([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],
                                [unit_vector1[1],unit_vector2[1],unit_vector3[1]],
                                [unit_vector1[2],unit_vector2[2],unit_vector3[2]]])
# [Add a description of what these quantities are and how they were obtained ... ]
test1_obs_time = [58577.489970740738,58583.545550740739,58590.545550740739]
sqr_mue = 1.7202098949957226E-002
tau_test1 = [0.22458337900908587,-0.10416868635938527,0.12041469264970059]

# [Add a description of what these position_vector* quantities are and how they were obtained ... ]
position_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
position_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
position_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
position_matrix1 = np.matrix([[position_vector1[0],position_vector2[0],position_vector3[0]],
                    [position_vector1[1],position_vector2[1],position_vector3[1]],
                    [position_vector1[2],position_vector2[2],position_vector3[2]]])

# [Add a description of what the various quantities below are and how they were obtained ... ]
inverse_matrix = [[-181.39997132286484,-49.141453222161964,732.7444648345615],
                [346.8208661801217,99.18166675238325,-1365.948759032001],
                [-167.22601164748204,-49.58309521064877,637.8813360694052]]
test1_dot_b =[[-0.005851587221372763],
[-0.001999529945903465],
[-0.0008667750934766205]]
test1_dot_a = [[0.005850249367257254],
[0.002025625186994423],
[0.0008769501443510741]]
test1_d_1 =1.0320244737743007
test1_d_2 = -1.0437988981068052
test1_r2pos = 1.003625082318933 
test1_dotr2 = 0.31100143213162895
test1_a1 = 0.5361692088746649
test1_b1 = 0.0032114744736032957
test1_a3 = 0.46383079112533504
test1_b3 = 0.0030602457022409226
test1_root2_r2m3 = 0.8528749073206177 # root2 1.054479693936294
test_roots1 = 1.245749369311985

def test_get_tau_values():
    """A test for the function `fortran_code_replica.get_tau_values` function
    
    Provided values
        test1_obs_time = [58577.489970740738,58583.545550740739,58590.545550740739]
        sqr_mue = 1.7202098949957226E-002
        
    Results:
        The code result are ->
        tau13 = 0.22458337900908587
        tau1 = -0.10416868635938527
        tau3 = 0.12041469264970059

        The fortran results are ->
        tau13:   0.22458337900908587 
        tau1:  -0.10416868635938527     
        tau3:   0.12041469264970059  
        
    """

    # [ I'd recommend restructuring the test a little to more clearly describe where the numbers come from ... ]
    # [ My text below is only a suggestion and you should read it to ensure that I am correct! ]
    # Expected values for tau13, tau1, tau3
    # These were calculated using the orbfit fortran code, using `test1_obs_time` as the input data
    # The orbfit fortran values for tau13, tau1, tau3 were printed out and then stored in the variable `tau_test1`
    expected_result = tau_test1

    # Call the function that we want to test
    result = fcr.get_tau_values(test1_obs_time,sqr_mue)

    # Check that the result matches the expected values
    assert np.allclose(result, expected_result, rtol=1e-9, atol=1e-9)
    # assert result1[0] == 0.22458337900908587
    # assert result1[1] == -0.10416868635938527
    # assert result1[2] == 0.12041469264970059

def test_make_a_b():
    """A test for the function make_a_b on the file fortran_code_replica.py, that provided tau. tau1, 
    tau3. test are performed for the return value of the function. 
    
    Provided values
        tau_test1 = [0.22458337900908587,-0.10416868635938527,0.12041469264970059]
    Results:
        The code result are ->
        a1  = 0.5361692088746649
        b1 = 0.0032114744736032957 
        a3 = 0.46383079112533504
        b3 = 0.0030602457022409226
        
        The fortran results are ->
        a1:   0.53616920887466490 
        b1:    3.2114744736032957E-003    
        a3:   0.46383079112533504 
        b3:    3.0602457022409226E-003
    """
    result1 = fcr.make_a_b(tau_test1[0],tau_test1[1],tau_test1[2])
    assert result1[0] == 0.53616920887466490 
    assert result1[1] == 3.2114744736032957E-003    
    assert result1[2] == 0.46383079112533504 
    assert result1[3] == 3.0602457022409226E-003
    
def test_get_d_and_pos_variables():
    """A test for the function get_d_and_pos_variables on the file fortran_code_replica.py, that 
    provided inverse matrix of the unit vector,rb(dot_b),ra(dot_a),the observer position
    matrix and the unit vector matrix.
    _ test are performed return value by the function. 
    
    Provided values
        inverse_matrix = [[-181.39997132286484,-49.141453222161964,732.7444648345615],
                [346.8208661801217,99.18166675238325,-1365.948759032001],
                [-167.22601164748204,-49.58309521064877,637.8813360694052]]
        test1_dot_b = [[-0.005851587221372763],
        [-0.001999529945903465],
        [-0.0008667750934766205]]
        test1_dot_a = [[0.005850249367257254],
        [0.002025625186994423],
        [0.0008769501443510741]] 
        postion_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
        postion_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
        postion_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
        position_matrix1 = [[postion_vector1[0],postion_vector2[0],postion_vector3[0]],
        [postion_vector1[1],postion_vector2[1],postion_vector3[1]],[postion_vector1[2],
        postion_vector2[2],postion_vector3[2]]]
        unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
        unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
        unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
        unit_vector_matrix1 = np.matrix([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],
        [unit_vector1[1],unit_vector2[1],unit_vector3[1]],[unit_vector1[2],
        unit_vector2[2],unit_vector3[2]]])
        
    Results:
        The code result are ->
        a2star = 1.0320244737743007
        b2star = -1.0437988981068052
        r22  = 1.003625082318933
        s2r2 = 0.31100143213162895
        
        The fortran results are ->
        a2star:    1.0320244737743007    
        b2star:   -1.0437988981068052     
        r22:    1.0036250823189330     
        s2r2:   0.31100143213162895
    """
    result1 = fcr.get_d_and_pos_variables(np.matrix(inverse_matrix),
                                        np.matrix(test1_dot_b),np.matrix(test1_dot_a),
                                        np.matrix(position_matrix1),np.matrix(unit_vector_matrix1))
    assert result1[0] == 1.0320244737743007
    assert result1[1] == -1.0437988981068052
    assert result1[2] == 1.0036250823189330 
    assert result1[3] == 0.31100143213162895

def test_generate_Cs():
    """A test for the function generate_Cs on the file fortran_code_replica.py, that provided with a2star(d1),b2star(d2),
    r22(r2pos),s2r2(dotr2). test are performed to check the return value of the function. 
    
    Provided values
        test1_d_1 =1.0320244737743007
        test1_d_2 = -1.0437988981068052
        test1_r2pos = 1.003625082318933 
        test1_dotr2 = 0.31100143213162895
        
    Results:
        The code result are ->
        c = -1.0895161396889805
        c3 = 2.8036979214270064
        c6 = -2.7106217754654516
        c8 = 1
        
        The fortran results are ->
        c: -1.0895161396889805     
        c3:  2.8036979214270064     
        c6: -2.7106217754654516     
        c8:  1.0000000000000000 
    """
    result1 = fcr.generate_Cs(test1_d_1,test1_d_2,test1_r2pos,test1_dotr2)

    assert len(result1) == 4
    assert result1[0] == 1.0000000000000000 
    assert result1[1] == -2.7106217754654516
    assert result1[2] == 2.8036979214270064
    assert result1[3] == -1.0895161396889805

def test_calculate_rhos():
    """A test for the function calculate_rhos on the file fortran_code_replica.py, that provided a1,b1,a3,b3,r2m3, 
    observer position matrix and the inverse of the unit vector matrix. test are performed for the return value of the function. 
    
    Provided values
        test1_a1 = 0.5361692088746649
        test1_b1 = 0.0032114744736032957
        test1_a3 = 0.46383079112533504
        test1_b3 = 0.0030602457022409226
        test1_r2m3 = 0.8528749073206177
        postion_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
        postion_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
        postion_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
        test1_pos_matrix = [[postion_vector1[0],postion_vector2[0],postion_vector3[0]],
        [postion_vector1[1],postion_vector2[1],postion_vector3[1]],[postion_vector1[2],
        postion_vector2[2],postion_vector3[2]]]
        inverse_matrix = [[-181.39997132286484,-49.141453222161964,732.7444648345615],
                [346.8208661801217,99.18166675238325,-1365.948759032001],
                [-167.22601164748204,-49.58309521064877,637.8813360694052]]
        (values come from test sample one root 1.054479693936294)
        
    Results:
        The code result are ->
        root = 1.054479693936294
        c1 = 0.5389081948687019
        c3 = 0.4664407978950121
        r2m3 = 0.8528749073206177
        gcap = 0.0008595774581503512     0.0003202762696971895    0.00013769941683439146
        crhom = -0.07076768208740036      0.14179458529004793     -0.07178811081637124
        rho1 = 0.13131676742203188
        rho3 = 0.1539061573094417 
        
        The fortran results are ->
        root:    1.0544796939362728     
        c_1:   0.53890819486870201     
        c_3:   0.46644079789501225     
        r2m3:   0.85287490732066884     
        gcap:    8.5957745815007369E-004   3.2027626969707845E-004   1.3769941683434983E-004
        crhom:   -7.0767682087375050E-002  0.14179458528999755       -7.1788110816345896E-002
        rho1:   0.13131676742198489     
        rho3:   0.15390615730938734 
        
    """
    result1 = fcr.calculate_rhos(test1_a1,test1_b1,test1_a3,
                                test1_b3,test1_root2_r2m3,position_matrix1, 
                                inverse_matrix)
    assert round(result1[0],13) == round(0.13131676742198489,13)
    assert round(result1[1],13) == round(0.15390615730938734,13)
    assert round(result1[2],13) == round(0.53890819486870201,13)
    assert round(result1[3],13) == round(0.46644079789501225,13)
    assert round(result1[4][0,0],14) == round(8.5957745815007369E-004,14)
    assert round(result1[4][1,0],14) == round(3.2027626969707845E-004,14)
    assert round(result1[4][2,0],14) == round(1.3769941683434983E-004,14)
    assert round(result1[5][0,0],13) == round(-7.0767682087375050E-002,13)
    assert round(result1[5][1,0],13) == round(0.14179458528999755,13)
    assert round(result1[5][2,0],12) == round(-7.1788110816345896E-002,12)
    assert round(result1[6],12) == round(0.85287490732066884,12)

def test_find_rho_and_r():
    """A test for the function find_rho_and_r on the file fortran_code_replica.py, that provided 
    the roots, a1, b1, a3, b3, the observation matrix, the inverse of the unit vector matrix, 
    and the unit vector matrix. test are performed on the return values. 
    
    Provided values
        test_roots1 = 1.245749369311985
        test_a1 = 0.5361692088746649
        test_b1 = 0.0032114744736032957
        test_a3 = 0.46383079112533504
        test_b3 = 0.0030602457022409226
        postion_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
        postion_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
        postion_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
        test_pos_matrix = [[postion_vector1[0],postion_vector2[0],postion_vector3[0]],
        [postion_vector1[1],postion_vector2[1],postion_vector3[1]],
        [postion_vector1[2],postion_vector2[2],postion_vector3[2]]]
        inverse_matrix = [[-181.39997132286484,-49.141453222161964,732.7444648345615],
                [346.8208661801217,99.18166675238325,-1365.948759032001],
                [-167.22601164748204,-49.58309521064877,637.8813360694052]]
        unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
        unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
        unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
        unit_vector_matrix1 = np.matrix([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],
        [unit_vector1[1],unit_vector2[1],unit_vector3[1]],[unit_vector1[2],unit_vector2[2],
        unit_vector3[2]]])
        
    Results:
        The code result are ->
        rho1 = 0.45894799661611524
        rho3 = 0.5326721981438207
        r1_1 = -1.211028021610874
        r1_2 = 0.16448007923237593
        r1_3 = -0.1303439561397366
        r3_1 = -1.242925080590317
        r3_2-0.019296651670140053
        r3_3-0.24064389301967798
        
        The fortran results are ->     
        rho1:   0.45894799661611524     
        rho3:   0.53267219814382072     
        r1_1:   -1.2110280216108740     
        r1_2:   0.16448007923237593     
        r1_3:  -0.13034395613973659     
        r3_1:   -1.2429250805903169     
        r3_2:   -1.9296651670140053E-002
        r3_3:  -0.24064389301967798   
        (a different root: 1.2457493693119874)
        
    """
    result1 = fcr.find_rho_and_r(test_roots1,test1_a1,
                                test1_b1,test1_a3,test1_b3,position_matrix1,
                                inverse_matrix,unit_vector_matrix1)    #return(rhos,r)
    assert len(result1) == 2
    assert type(result1[0]) and type(result1[1])== list
    assert result1[0][0] == 0.45894799661611524 
    assert result1[0][1] == 0.53267219814382072
    assert result1[1][0][0,0] == -1.2110280216108740
    assert result1[1][0][1,0] == 0.16448007923237593 
    assert result1[1][0][2,0] == -0.13034395613973659
    assert result1[1][1][0,0] == -1.2429250805903169 
    assert result1[1][1][1,0] == -1.9296651670140053E-002
    assert result1[1][1][2,0] == -0.24064389301967798

def test_gauss_method_8th():
    """A test for the function gauss_method_8th on the file fortran_code_replica.py, that provided the 
    observer position matrix, the unit vector matrix, and the observation time.
    a tets to check the returned values of the function is performed
    
    Provided values
        postion_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
        postion_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
        postion_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
        test_pos_matrix  = [[postion_vector1[0],postion_vector2[0],postion_vector3[0]],
        [postion_vector1[1],postion_vector2[1],postion_vector3[1]],[postion_vector1[2],
        postion_vector2[2],postion_vector3[2]]] From test sample 1
        test1_obs_time = test1_obs_time = [58577.489970740738,
        58583.545550740739,58590.545550740739] From test sample1 
        unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
        unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
        unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
        test1_unit_vector = np.matrix([[unit_vector1[0],unit_vector2[0],
        unit_vector3[0]],[unit_vector1[1],unit_vector2[1],unit_vector3[1]],
        [unit_vector1[2],unit_vector2[2],unit_vector3[2]]]) From test sample 1

    Results:
        The code result are ->
        root1 = 1.245749369311985
        root2 = 1.054479693936294
        root3 = 0.9836605218148519
        rho values
        root1 rho1 = 0.45894799661611524
        root1 rho3 = 0.5326721981438207
        root2 rho1 = 0.13131676742203188
        root2 rho3 = 0.1539061573094417
        root3 rho1 = -0.061154610273003375
        root3 rho3 = -0.06853359800812431
        r values
        r1
        root1 r1_1 = -1.211028021610874
        root1 r1_2 = 0.16448007923237593
        root1 r1_3 = -0.1303439561397366
        root2 r1_1 = -1.0387491242844016
        root2 r1_2 = -0.11319992500130222
        root2 r1_3 = -0.10676389829448926
        root3 r1_1 = -0.9375415660251583
        root3 r1_2 = -0.27632678291617774
        root3 r1_3 = -0.09291147576820177
        r3
        root1 r3_1 = -1.242925080590317
        root1 r3_2 = -0.019296651670140053
        root1 r3_3 = -0.24064389301967798
        root2 r3_1 = -0.995179394152981
        root2 r3_2 = -0.3026587102533421
        root2 r3_3 = -0.19831503411073065
        root3 r3_1 = -0.8496845927956426
        root3 r3_2 = -0.46907010985685227
        root3 r3_3 = -0.17345636146458615
        
        The fortran results are ->
        root1:    1.2457493693119874 
        root2:    1.0544796939362728  
        root3:   0.98366052181486430  
        rho values
        root1 rho1:   0.45894799661611524     
        root1 rho3:   0.53267219814382072 
        root2 rho1:   0.13131676742198489     
        root2 rho3:   0.15390615730938734   
        root3 rho1:   -6.1154610272977687E-002
        root3 rho3:   -6.8533598008092847E-002
        r values
        r1
        root1 r1_1:   -1.2110280216108740     
        root1 r1_2:   0.16448007923237593     
        root1 r1_3:  -0.13034395613973659  
        root2 r1_1:   -1.0387491242843767     
        root2 r1_2:  -0.11319992500134204     
        root2 r1_3:  -0.10676389829448588   
        root3 r1_1:  -0.93754156602517169     
        root3 r1_2:  -0.27632678291615598     
        root3 r1_3:   -9.2911475768203619E-002
        r3
        root1 r3_1:  -1.2429250805903169     
        root1 r3_2:  -1.9296651670140053E-002
        root1 r3_3:  -0.24064389301967798 
        root2 r3_1:  -0.99517939415294532     
        root2 r3_2:  -0.30265871025338281     
        root2 r3_3:  -0.19831503411072457   
        root3 r3_1:  -0.84968459279566322     
        root3 r3_2:  -0.46907010985682873     
        root3 r3_3:  -0.17345636146458968   
    """
    results1 = fcr.gauss_method_8th(np.matrix(position_matrix1),
                                    test1_obs_time,np.matrix(unit_vector_matrix1))
    assert len(results1) == 3
    assert len(results1[0]) == 3
    #roots
    assert round(results1[0][0],13) == round( 1.2457493693119874 ,13)
    assert round(results1[1][0],13) == round(1.0544796939362728 ,13)
    assert round(results1[2][0],13) == round(0.98366052181486430,13)
    
    # rhos per root
    assert results1[0][1][0] == 0.45894799661611524
    assert results1[0][1][1] == 0.53267219814382072
    
    assert round(results1[1][1][0],13) == round(0.13131676742198489,13)
    assert round(results1[1][1][1],13) == round(0.15390615730938734,13)
    
    assert round(results1[2][1][0],13) == round(-6.1154610272977687E-002,13)
    assert round(results1[2][1][1],13) == round(-6.8533598008092847E-002,13)
    
    # r values
    # r_1 
    assert results1[0][2][0][0,0] == -1.2110280216108740 
    assert results1[0][2][0][1,0] == 0.16448007923237593 
    assert results1[0][2][0][2,0] == -0.13034395613973659
    
    assert round(results1[1][2][0][0,0],13) == round(-1.0387491242843767 ,13)
    assert round(results1[1][2][0][1,0],13) == round(-0.11319992500134204,13)
    assert round(results1[1][2][0][2,0],14) == round(-0.10676389829448588,14)
    
    assert round(results1[2][2][0][0,0],13) == round(-0.93754156602517169,13)
    assert round(results1[2][2][0][1,0],13) == round(-0.27632678291615598 ,13)
    assert round(results1[2][2][0][2,0],14) == round(-9.2911475768203619E-002,14)
    
    # r_3 values
    assert results1[0][2][1][0,0] ==-1.2429250805903169 
    assert results1[0][2][1][1,0] ==-1.9296651670140053E-002 
    assert results1[0][2][1][2,0] ==-0.24064389301967798 
    
    assert round(results1[1][2][1][0,0],12) == round(-0.99517939415294532,12)
    assert round(results1[1][2][1][1,0],12) == round(-0.30265871025338281 ,12)
    assert round(results1[1][2][1][2,0],13) == round(-0.19831503411072457,13)
    
    assert round(results1[2][2][1][0,0],12) == round(-0.84968459279566322 ,12)
    assert round(results1[2][2][1][1,0],12) == round(-0.46907010985682873 ,12)
    assert round(results1[2][2][1][2,0],14) == round(-0.17345636146458968,14)
