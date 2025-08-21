import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../Honor_thesis_Code")))
from gauss import gauss_method as gm
np.set_printoptions(threshold=np.inf,precision=20)

test1_ra_dec = [2.1260971174823160,-7.2033616739254860E-002]
test2_ra_dec =[2.2018577151544165,-9.1656450482163324E-002]
test3_ra_dec = [2.2892342062190609,-0.11198856664053504]
test1_times = [58577.489970740738,58583.545550740739,58590.545550740739]
test_d_t1 = -6.055580000000191
test_d_t3 = 7
testd_t = 13.055580000000191
mu = 1.7202098949957226E-002
unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
test1_rho = np.array([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],[unit_vector1[1],unit_vector2[1],unit_vector3[1]],[unit_vector1[2],unit_vector2[2],unit_vector3[2]]])
test1_values_ra_dec = [test1_ra_dec,test2_ra_dec,test3_ra_dec]
test1_a_1 = 0.5361692088746649
test1_b_1 = 10.85279479419046
test1_a_3 = 0.4638307911253351
test1_b_3 = 10.341735205810208
postion_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
postion_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
postion_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
test1_position_R = np.array([[postion_vector1[0],postion_vector2[0],postion_vector3[0]],[postion_vector1[1],postion_vector2[1],postion_vector3[1]],[postion_vector1[2],postion_vector2[2],postion_vector3[2]]])
test1_B_matrix = np.array([[115.62989788917604,85.74605554189431,50.08435887453609],
[-225.65321388125662,-170.38814074716754,-104.27900008385723],
[ 111.21597772750607,85.55454550094318,54.77127152627216]])
test1_d_1 =1.0320244737742499 
test1_d_2 = -3527.3938312967603
test1_magnitude_R_2 = 1.0018109014773862 
test1_pos_2_dot_p2 = 0.31100143213162895
test1_1_root =1.0544796939365233

def test_calculate_unit_vector():
    """
    A test for the function calculate_unit_vector on the file gauss_method.py, that provided the Right Ascencion (RA) and 
    Declination (DEC) observation, provieds the Equritorial helocentric postion of the object.
    Three test are performend each with a different RA and DEC values, to check the return 
    position coordinate by the function. 
    
    Provided values
        test1_ra_dec = [2.1260971174823160,-7.2033616739254860E-002]
        test2_ra_dec =[2.2018577151544165,-9.1656450482163324E-002]
        test3_ra_dec = [2.2892342062190609,-0.11198856664053504]
        First value on the provided array is the RA, and the second value is the DEC
        
    Results:
        The code result are ->
        [-0.5258317338986549, 0.8475382670837678, -0.07197133772397166]
        [-0.5875255158277471, 0.8040126628785391, -0.09152817152276324]
        [-0.6540863217081979, 0.7481189653618031, -0.11175463041961636]
        
        The fortran results are ->
        vector 1 [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
        vector 2 [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
        vector 3 [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
        
    Due to numpy matrice data type claulation are only precide to the 17 digit. 
    """
    result1 = gm.calculate_unit_vector(test1_ra_dec[0],test1_ra_dec[1])
    result2 = gm.calculate_unit_vector(test2_ra_dec[0],test2_ra_dec[1])
    result3 = gm.calculate_unit_vector(test3_ra_dec[0],test3_ra_dec[1])
    result4 = gm.calculate_unit_vector([test1_ra_dec[0],test2_ra_dec[0],test3_ra_dec[0]],[test1_ra_dec[1],test2_ra_dec[1],test2_ra_dec[1]])
    assert len(result1) == len(result2) == len(result3) == 3
    assert result1.all() == np.array([-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]).all()
    assert result2.all() == np.array([-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]).all()
    assert result3.all() == np.array([-0.65408632170819792,0.74811896536180311,-0.11175463041961636]).all()
    assert result4.all() == np.array(([[-0.5258317338986549, 0.8475382670837678, -0.07197133772397166],
                                    [-0.5875255158277471, 0.8040126628785391, -0.09152817152276324],
                                    [-0.6540863217081979, 0.7481189653618031, -0.11175463041961636]])).all()
test_calculate_unit_vector()

def test_calculate_delta_t():
    """A test for the function calculate_delta_t on the file gauss_method.py, that provided the time time 
    of observations, as a list of three length, an provides the the three tau values used on the gauss method.
    
    Provided Values:
        test1_times = [58577.489970740738,58583.545550740739,58590.545550740739]
    Results:
        The code result are ->
            tau13 = 13.055580000000191      times mue = 0.22458337900908584
            tau1 = -6.055580000000191      times mue = -0.10416868635938527
            tau3 = 7.0                     times mue = 0.12041469264970059
        
        The fortran results are ->
            tau13 =  0.22458337900908587 
            tau1 =  -0.10416868635938527     
            tau3 =   0.12041469264970059  
    """
    result = gm.calculate_delta_t(test1_times)
    # Multiply by 1.7202098949957226E-002 (mue square), to compare with the fortran values
    assert result[0] * 1.7202098949957226E-002  == -0.10416868635938527,15 # -6.055580000000191
    assert result[1] * 1.7202098949957226E-002 == 0.12041469264970059  # 7
    assert round(result[2] * 1.7202098949957226E-002,15) == round(0.22458337900908587,15) # 13.055580000000191

def test_calculate_a_and_b():
    """A test for the function calculate_a_and_b on the file gauss_method.py, the tau values of the gauss method it
    generates the nessecary a_1, b_1, a_3, and b_3 values, this test is to make the function works correctly 
    and is capable of calculate the expected values
    
    Provided Values:
        test_d_t1 = -6.055580000000191
        test_d_t3 = 7
        testd_t = 13.055580000000191
        mu = 1.7202098949957226E-002
        
    Results:
        The code result are ->
        a1 = 0.5361692088746649  
        b1 = 10.85279479419046    times mu square = 0.0032114744736032952
        a3 = 0.4638307911253351  
        b3 = 10.341735205810208   times mu square = 0.0030602457022409217
        
        The fortran results are ->
            Edit a1:   0.53616920887466490 
            Edit b1:    3.2114744736032957E-003    
            Edit a3:   0.46383079112533504 
            Edit b3:    3.0602457022409226E-003
    """
    result = gm.calculate_a_and_b(test_d_t1,test_d_t3,testd_t)
    
    # Some values are multiplied by mu square (result [1] and [3]) to match the scaled values from the Fortran code.
    assert result[0] == 0.53616920887466490 
    assert round(result[1] * mu** 2,17) == round(3.2114744736032957E-003,17)
    assert round(result[2],15) == round(0.46383079112533504,15)
    assert round(result[3] *mu **2,17) == round(3.0602457022409226E-003,17)
    
def test_get_unit_matrix():
    """A test for the function get_unit_matrix on the file gauss_method.py, a list with 3 sets of RA and DEC values
    provides the unit matrix for the observation.
    A matrix is made with the fortran results named test1_rho and used to compare the results 
    from the function to the fortran code. 
    
    Provided Values:
        test1_ra_dec = [2.1260971174823160,-7.2033616739254860E-002]
        test2_ra_dec =[2.2018577151544165,-9.1656450482163324E-002]
        test3_ra_dec = [2.2892342062190609,-0.11198856664053504]
        test1_values_ra_dec = [test1_ra_dec,test2_ra_dec,test3_ra_dec]
        
    Results:
        The code result are ->
        [[-0.5258317338986549  -0.5875255158277471  -0.6540863217081979 ]
        [ 0.8475382670837678   0.8040126628785391   0.7481189653618031 ]
        [-0.07197133772397166 -0.09152817152276324 -0.11175463041961636]]
        
        The fortran results are ->
        [[-0.5258317338986549  -0.5875255158277471  -0.6540863217081979 ]
        [ 0.8475382670837678   0.8040126628785391   0.7481189653618031 ]
        [-0.07197133772397166 -0.09152817152276324 -0.11175463041961636]]
            
    """
    result1 = gm.get_unit_matrix(test1_values_ra_dec)
    assert type(result1) == np.ndarray 
    assert (result1 == test1_rho).all
    print(result1)

# Test for function calculate_postion are on file test_wis

def test_calculate_values_for_cs():
    """A test for the function calculate_values_for_cs on the file gauss_method.py, This function provided it'sargument; 
    the a1, b1, a3, b3, position matrix R, unit matrix of the object, and the "B" matrix 
    provides values neccesary to calculate the four coefficient used on the eight degreee equation used on the gauss method
    Provided Values:
        test1_a_1 = 0.5361692088746649
        test1_b_1 = 10.85279479419046
        test1_a_3 = 0.4638307911253351
        test1_b_3 = 10.341735205810208
        postion_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
        postion_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
        postion_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
        test1_position_R = np.matrix([[postion_vector1[0],postion_vector2[0],postion_vector3[0]],[postion_vector1[1],postion_vector2[1],postion_vector3[1]],[postion_vector1[2],postion_vector2[2],postion_vector3[2]]])
        unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
        unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
        unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
        test1_rho = np.matrix([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],[unit_vector1[1],unit_vector2[1],unit_vector3[1]],[unit_vector1[2],unit_vector2[2],unit_vector3[2]]])
        test1_B_matrix = np.matrix([[115.62989788917604,85.74605554189431,50.08435887453609],
        [-225.65321388125662,-170.38814074716754,-104.27900008385723],
        [ 111.21597772750607,85.55454550094318,54.77127152627216]])
        
    Results:
        The code result are ->
        d1 = 1.0320244737742499 
        d_2 = -3527.3938312967603 
        Magnitude r2 = 1.0018109014773862 
        dot product = 0.31100143213162895
        
        The fortran results are ->
        No Fortran result is provided, because this function performs mathematical operations differently 
        than the Fortran orbfit code.  
            
    """
    result = gm.calculate_values_for_cs(test1_a_1,test1_b_1,test1_a_3,test1_b_3,test1_position_R,test1_rho,test1_B_matrix)
    assert result[0] == 1.0320244737742499
    assert result[1] == -3527.3938312967603
    assert result[2] == 1.0018109014773862
    assert result[3] == 0.31100143213162895

def test_calculate_cs():
    """A test for the function calculate_cs on the file gauss_method.py, with d1, d2, magnitude of R2, and the 
    dot product it calculates the four coefficient used on the eight degreee equation for the gauss method
    Provided Values:
        test1_d_1 =1.0320244737742499 
        test1_d_2 = -3527.3938312967603
        test1_magnitude_R_2 = 1.0018109014773862 
        test1_pos_2_dot_p2 = 0.31100143213162895
        
    Results:
        The code result are ->
        c = 1.0895161396889808
        c3 =-2.8036979214269007
        c6 = 2.7106217754653157
        c8 = -1
        
        The fortran results are ->
        c:   -1.0895161396889805     
        c3:    2.8036979214270064     
        c6:   -2.7106217754654516     
        c8:    1.0000000000000000 
        
    """
    result = gm.calculate_cs(test1_d_1,test1_d_2,test1_magnitude_R_2,test1_pos_2_dot_p2)
    # Assertion results are multiplied by -1 to match the sign convention used in the Fortran code.
    assert round(-1*result[0],14) == round(-1.0895161396889805,14)
    assert round(-1*result[1],12) == round(2.8036979214270064,12)
    assert round(-1*result[2],12) == round(-2.7106217754654516,12)
    assert round(-1*result[3],14) == 1.0000000000000000

def test_calculate_rhos_and_r():
    """A test for the function calculate_rhos_and_r on the file gauss_method.py, with d1, d2, magnitude of R2, and the 
    dot product it calculates the four coefficient used on the eight degreee equation for the gauss method
    Provided Values:
        test1_1_sigma = 0.00025237609721519305 # for root 1.0544796939365233
        test1_B_matrix = np.matrix([[115.62989788917604,85.74605554189431,50.08435887453609],
        [-225.65321388125662,-170.38814074716754,-104.27900008385723],
        [ 111.21597772750607,85.55454550094318,54.77127152627216]])
        test1_a_1 = 0.5361692088746649
        test1_b_1 = 10.85279479419046
        test1_a_3 = 0.4638307911253351
        test1_b_3 = 10.341735205810208
        postion_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
        postion_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
        postion_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
        test1_position_R = np.matrix([[postion_vector1[0],postion_vector2[0],postion_vector3[0]],[postion_vector1[1],postion_vector2[1],postion_vector3[1]],[postion_vector1[2],postion_vector2[2],postion_vector3[2]]])
        unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
        unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
        unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
        test1_rho = np.matrix([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],[unit_vector1[1],unit_vector2[1],unit_vector3[1]],[unit_vector1[2],unit_vector2[2],unit_vector3[2]]])
        
    Results:
        The code result are ->
        root = 1.0544796939365233
        # p1 0.13131676742254747
        # p2 0.1417945852906265
        # p3 0.15390615731009605
        # r11-1.0387491242846727
        # r12 -0.11319992500086523
        # r13 -0.10676389829452637
        # r21 -1.0239826861472827
        # r22 -0.20217673754170323
        # r23 -0.1500381624522676
        # r31 -0.995179394153409
        # r32 -0.3026587102528526
        # r33 -0.19831503411080378
        
        The fortran results are ->
        root:    1.0544796939362728
        p1:   0.13131676742198489
        p2:   0.14179458528999755
        p3:   0.15390615730938734
        r11:   -1.0387491242843767
        r12:  -0.11319992500134204
        r13:  -0.10676389829448588
        r21:   -1.0239826861469132
        r22:  -0.20217673754220894
        r23:  -0.15003816245221002
        r31:  -0.99517939415294532
        r32:  -0.30265871025338281
        r33:  -0.19831503411072457

    sigma provided value depends on the given root
    """
    result1 = gm.calculate_rhos_and_r(test1_1_root,test1_B_matrix,test1_a_1,test1_b_1,test1_a_3,test1_b_3,test1_position_R,test1_rho)
    assert len(result1) == 3 and len(result1[1]) == 3 and type(result1[2][0][0]) == np.float64
    assert result1[0] == test1_1_root # root
    assert round(result1[1][0],11)  == round(0.13131676742198489,11) # rho 1
    assert round(result1[1][1],11)  == round(0.14179458528999755,11)# rho 2
    assert round(result1[1][2],11)  ==  round(0.15390615730938734,11)# rho 3
    assert round(result1[2][0,0],11)  == round(-1.0387491242843767,11)# matrix component r_1 (1,0)
    assert round(result1[2][1,0],12)  == round(-0.11319992500134204,12)# matrix component r_1 (2,0)
    assert round(result1[2][2,0],11)  == round( -0.10676389829448588,11)# matrix component r_1 (3,0)
    assert round(result1[2][0,1],12)  == round(-1.0239826861469132,12)# matrix component r_2 (1,0)
    assert round(result1[2][1,1],12)  == round(-0.20217673754220894,12)# matrix component r_2 (2,0)
    assert round(result1[2][2,1],12)  == round(-0.15003816245221002,12)# matrix component r_2 (3,0)
    assert round(result1[2][0,2],12)  == round(-0.99517939415294532,12)# matrix component r_3 (1,0)
    assert round(result1[2][1,2],12)  == round( -0.30265871025338281,12)# matrix component r_3 (2,0)
    assert round(result1[2][2,2],12) == round(-0.19831503411072457,12)# matrix component r_2 (3,0)    

def test_gauss_method():
    """A test for the function  on the file gauss_method.py, this function provided three observation times, the rho matrix and psotion
    of the observer returns the rho and r values using the gauss method (eight degree equation)
    Provided Values:
        test1_times = [58577.489970740738,58583.545550740739,58590.545550740739]
        unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
        unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
        unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
        test1_rho = np.matrix([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],[unit_vector1[1],unit_vector2[1],unit_vector3[1]],[unit_vector1[2],unit_vector2[2],unit_vector3[2]]])
        postion_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
        postion_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
        postion_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
        test1_position_R = np.matrix([[postion_vector1[0],postion_vector2[0],postion_vector3[0]],[postion_vector1[1],postion_vector2[1],postion_vector3[1]],[postion_vector1[2],postion_vector2[2],postion_vector3[2]]])

    Results:
        The code result are ->
        root = 0.9836605218147869
        p1 = -0.06115461027327507
        p2 = -0.06465851941455836
        p3 = -0.06853359800837579
        r11 = -0.9375415660250154
        r12 = -0.276326782916408
        r13 = -0.09291147576818222
        r21 = -0.9026862193111291
        r22 = -0.36816764801526075
        r23 = -0.13114188727340442
        r31 = -0.8496845927954781
        r32 = -0.4690701098570404
        r33 = -0.17345636146455806

        root = 1.0544796939365233
        p1 = 0.13131676742254747
        p2 = 0.1417945852906265
        p3 = 0.15390615731009605
        r11 = -1.0387491242846727
        r12 = -0.11319992500086523
        r13 = -0.10676389829452637
        r21 = -1.0239826861472827
        r22 = -0.20217673754170323
        r23 = -0.1500381624522676
        r31 = -0.995179394153409
        r32 = -0.3026587102528526
        r33 = -0.19831503411080378

        root = 1.2457493693117867
        p1 = 0.4589479966158301
        p2 = 0.4921102037009649
        p3 = 0.5326721981435508
        r11 = -1.211028021610724
        r12 = 0.16448007923213429
        r13 = -0.13034395613971606
        r21 = -1.2298020505563332
        r22 = 0.07948145566433512
        r23 = -0.18210191046123192
        r31 = -1.2429250805901404
        r32-0.019296651670341947
        r33 = -0.24064389301964784
        
        The fortran results are ->
        root:   0.98366052181486430
        p1:   -6.1154610272977687E-002
        p2:   -6.4658519414266299E-002
        p3:   -6.8533598008092847E-002
        r11:  -0.93754156602517169
        r12:  -0.27632678291615598
        r13:   -9.2911475768203619E-002
        r21:  -0.90268621931130066
        r22:  -0.36816764801502594
        r23:  -0.13114188727343115
        r31:  -0.84968459279566322
        r32:  -0.46907010985682873
        r33:  -0.17345636146458968

        root:    1.0544796939362728
        p1:   0.13131676742198489
        p2:   0.14179458528999755
        p3:   0.15390615730938734
        r11:   -1.0387491242843767
        r12:  -0.11319992500134204
        r13:  -0.10676389829448588
        r21:   -1.0239826861469132
        r22:  -0.20217673754220894
        r23:  -0.15003816245221002
        r31:  -0.99517939415294532
        r32:  -0.30265871025338281
        r33:  -0.19831503411072457

        root:    1.2457493693119874
        p1:   0.45894799661611524
        p2:   0.49211020370124259
        p3:   0.53267219814382072
        r11:   -1.2110280216108740
        r12:   0.16448007923237593
        r13:  -0.13034395613973659
        r21:   -1.2298020505564962
        r22:    7.9481455664558387E-002
        r23:  -0.18210191046125732
        r31:   -1.2429250805903169
        r32:   -1.9296651670140053E-002
        r33:  -0.24064389301967798
        
    """
    result = gm.gauss_method(test1_times,test1_rho,test1_position_R)
    assert len(result) == 3 
    assert round(result[0][0],11) == round(1.2457493693119874,11)# root
    assert round(result[0][1][0],11)  == round(0.45894799661611524,11) # rho 1
    assert round(result[0][1][1],11)  == round(0.49211020370124259,11)# rho 2
    assert round(result[0][1][2],11)  ==  round(0.53267219814382072,11)# rho 3
    assert round(result[0][2][0][0,0],11)  == round(-1.2110280216108740,11)# matrix component r_1 (1,0)
    assert round(result[0][2][0][1,0],12)  == round(0.16448007923237593,12)# matrix component r_1 (2,0)
    assert round(result[0][2][0][2,0],11)  == round(-0.13034395613973659,11)# matrix component r_1 (3,0)
    assert round(result[0][2][1][0,0],12)  == round(-1.2298020505564962,12)# matrix component r_2 (1,0)
    assert round(result[0][2][1][1,0],11)  == round(7.9481455664558387E-002,11)# matrix component r_2 (2,0)
    assert round(result[0][2][1][2,0],12)  == round(-0.18210191046125732,12)# matrix component r_2 (3,0)
    assert round(result[0][2][2][0,0],12)  == round(-1.2429250805903169,12)# matrix component r_3 (1,0)
    assert round(result[0][2][2][1,0],12)  == round(-1.9296651670140053E-002,12)# matrix component r_3 (2,0)
    assert round(result[0][2][2][2,0],12) == round(-0.24064389301967798,12)# matrix component r_2 (3,0)
    
    assert round(result[1][0],11) == round(1.0544796939362728,11)# root
    assert round(result[1][1][0],11)  == round(0.13131676742198489,11) # rho 1
    assert round(result[1][1][1],11)  == round(0.14179458528999755,11)# rho 2
    assert round(result[1][1][2],11)  ==  round(0.15390615730938734,11)# rho 3
    assert round(result[1][2][0][0,0],11)  == round(-1.0387491242843767,11)# matrix component r_1 (1,0)
    assert round(result[1][2][0][1,0],12)  == round(-0.11319992500134204,12)# matrix component r_1 (2,0)
    assert round(result[1][2][0][2,0],11)  == round( -0.10676389829448588,11)# matrix component r_1 (3,0)
    assert round(result[1][2][1][0,0],12)  == round(-1.0239826861469132,12)# matrix component r_2 (1,0)
    assert round(result[1][2][1][1,0],12)  == round(-0.20217673754220894,12)# matrix component r_2 (2,0)
    assert round(result[1][2][1][2,0],12)  == round(-0.15003816245221002,12)# matrix component r_2 (3,0)
    assert round(result[1][2][2][0,0],12)  == round(-0.99517939415294532,12)# matrix component r_3 (1,0)
    assert round(result[1][2][2][1,0],12)  == round( -0.30265871025338281,12)# matrix component r_3 (2,0)
    assert round(result[1][2][2][2,0],12) == round(-0.19831503411072457,12)# matrix component r_2 (3,0)
    
    assert round(result[2][0],11) == round(0.98366052181486430,11)# root
    assert round(result[2][1][0],11)  == round(-6.1154610272977687E-002,11) # rho 1
    assert round(result[2][1][1],11)  == round(-6.4658519414266299E-002,11)# rho 2
    assert round(result[2][1][2],11)  ==  round(-6.8533598008092847E-002,11)# rho 3
    assert round(result[2][2][0][0,0],11)  == round(-0.93754156602517169,11)# matrix component r_1 (1,0)
    assert round(result[2][2][0][1,0],12)  == round(-0.27632678291615598,12)# matrix component r_1 (2,0)
    assert round(result[2][2][0][2,0],11)  == round(-9.2911475768203619E-002,11)# matrix component r_1 (3,0)
    assert round(result[2][2][1][0,0],12)  == round(-0.90268621931130066,12)# matrix component r_2 (1,0)
    assert round(result[2][2][1][1,0],12)  == round(-0.36816764801502594,12)# matrix component r_2 (2,0)
    assert round(result[2][2][1][2,0],12)  == round(-0.13114188727343115,12)# matrix component r_2 (3,0)
    assert round(result[2][2][2][0,0],11)  == round(-0.84968459279566322,11)# matrix component r_3 (1,0)
    assert round(result[2][2][2][1,0],12)  == round(-0.46907010985682873,12)# matrix component r_3 (2,0)
    assert round(result[2][2][2][2,0],12) == round(-0.17345636146458968,12)# matrix component r_2 (3,0)

# test_gauss_method()
print("hello")