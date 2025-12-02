import sys
import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../Honor_thesis_Code")))
from gauss import gauss_method as gm
np.set_printoptions(threshold=np.inf,precision=20)

# Test values aqure from sample one observations and fortran
test1_ra_dec = [2.1260971174823160,-7.2033616739254860E-002]
test2_ra_dec =[2.2018577151544165,-9.1656450482163324E-002]
test3_ra_dec = [2.2892342062190609,-0.11198856664053504]
test1_values_ra_dec = np.array([[test1_ra_dec[0],test2_ra_dec[0],test3_ra_dec[0]],
                            [test1_ra_dec[1],test2_ra_dec[1],test3_ra_dec[1]]])
test1_times = [58577.489970740738,58583.545550740739,58590.545550740739]
mu = 1.7202098949957226E-002    # value from Fortran code
test1_rho =np.array([[-0.52583173389865490, -0.58752551582774715, -0.65408632170819792 ],   # x values
                    [ 0.84753826708376778 , 0.80401266287853912, 0.74811896536180311 ], # y values
                    [-7.1971337723971657E-002,-9.1528171522763241E-002, -0.11175463041961636]]) # Z values

# print(test1_rho.T)
test1_position_R = [[-0.96969860078090808, -0.94067474928282591 , -0.89451148183020490 ],
                            [-0.22449591050121329, -0.31618137964297799 ,  -0.41779882542249258 ],
                            [-9.7312854877537963E-002 ,-0.13705996332878803 ,-0.18111530838131065]]

# test input values from running python code that produced the expected result
test_d_t1 = -6.055580000000191
test_d_t3 = 7
testd_t = 13.055580000000191
test1_a_1 = 0.5361692088746649
test1_b_1 = 10.85279479419046
test1_a_3 = 0.4638307911253351
test1_b_3 = 10.341735205810208
test1_B_matrix = np.array([[115.62989788917604,85.74605554189431,50.08435887453609],
[-225.65321388125662,-170.38814074716754,-104.27900008385723],
[ 111.21597772750607,85.55454550094318,54.77127152627216]])
test1_d_1 =1.0320244737742499 
test1_d_2 = -3527.3938312967603
test1_magnitude_R_2 = 1.0018109014773862 
test1_pos_2_dot_p2 = 0.31100143213162895
test1_1_root =1.0544796939365233


# Expected Results:
# test_calculate_unit_vector
exp_unit_vec1 = np.array([[-0.52583173389865490],[0.84753826708376778],[-7.1971337723971657E-002]]) # Expected result of unit vector 1
exp_unit_vec2 = np.array([[-0.58752551582774715],[0.80401266287853912],[-9.1528171522763241E-002]])
exp_unit_vec3 = np.array([[-0.65408632170819792],[0.74811896536180311],[-0.11175463041961636]])
exp_unit_vec4 = np.array(([[-0.5258317338986549, -0.58752551582774715, -0.65408632170819792],
                                    [0.84753826708376778, 0.80401266287853912, 0.74811896536180311],
                                    [-7.1971337723971657E-002, -9.1528171522763241E-002, -0.11175463041961636]]))

# test_calculate_delta_t
exp_cal_dt1_0 = -0.10416868635938527
exp_cal_dt1_1 = 0.12041469264970059  
exp_cal_dt1_2 = 0.22458337900908587 # is rounded to 15

# test_calculate_a_and_b
exp_cal_ab1_0 = 0.53616920887466490 
exp_cal_ab1_1 = 3.2114744736032957E-003
exp_cal_ab1_2 = 0.46383079112533504
exp_cal_ab1_3 = 3.0602457022409226E-003

# test_calculate_values_for_cs
exp_cal_valcs1_0 = 1.0320244737742499
exp_cal_valcs1_1 = -3527.3938312967603
exp_cal_valcs1_2 = 1.0018109014773862
exp_cal_valcs1_3 = 0.31100143213162895

# test_calculate_rhos_and_r
exp_calc_rhosr1_10 = 0.13131676742198489
exp_calc_rhosr1_11 = 0.14179458528999755 # rho 2
exp_calc_rhosr1_12 =  0.15390615730938734 # rho 3
exp_calc_rhosr1_200 = -1.0387491242843767 # matrix component r_1 (1,0)
exp_calc_rhosr1_210 = -0.11319992500134204 # matrix component r_1 (2,0)
exp_calc_rhosr1_220 =  -0.10676389829448588 # matrix component r_1 (3,0)
exp_calc_rhosr1_201 = -1.0239826861469132 # matrix component r_2 (1,0)
exp_calc_rhosr1_211= -0.20217673754220894 # matrix component r_2 (2,0)
exp_calc_rhosr1_221= -0.15003816245221002 # matrix component r_2 (3,0)
exp_calc_rhosr1_202= -0.99517939415294532 # matrix component r_3 (1,0)
exp_calc_rhosr1_212= -0.30265871025338281 # matrix component r_3 (2,0) # need to fix
exp_calc_rhosr1_222= -0.19831503411072457 # matrix component r_2 (3,0)

# test_gauss_method
exp_gasuss1_00 = 1.2457493693119874 # root
exp_gasuss1_2_000 = -1.2110280216108740 # matrix component r_1 (1,0)
exp_gasuss1_2_010 = 0.16448007923237593 # matrix component r_1 (2,0)
exp_gasuss1_2_020 = -0.13034395613973659 # matrix component r_1 (3,0)
exp_gasuss1_2_001 = -1.2298020505564962 # matrix component r_2 (1,0)
exp_gasuss1_2_011 = 7.9481455664558387E-002 # matrix component r_2 (2,0)
exp_gasuss1_2_021 = -0.18210191046125732 # matrix component r_2 (3,0)
exp_gasuss1_2_002 = -1.2429250805903169 # matrix component r_3 (1,0)
exp_gasuss1_2_012 = -1.9296651670140053E-002 # matrix component r_3 (2,0)
exp_gasuss1_2_022 = -0.24064389301967798 # matrix component r_2 (3,0)
exp_gasuss1_01 = 1.0544796939362728 # root
exp_gasuss1_2_100 = -1.0387491242843767 # matrix component r_1 (1,0)
exp_gasuss1_2_110 = -0.11319992500134204 # matrix component r_1 (2,0)
exp_gasuss1_2_120 = -0.10676389829448588 # matrix component r_1 (3,0)
exp_gasuss1_2_101 = -1.0239826861469132 # matrix component r_2 (1,0)
exp_gasuss1_2_111 = -0.20217673754220894 # matrix component r_2 (2,0)
exp_gasuss1_2_121 = -0.15003816245221002 # matrix component r_2 (3,0)
exp_gasuss1_2_102 = -0.99517939415294532 # matrix component r_3 (1,0)
exp_gasuss1_2_112 = -0.30265871025338281# matrix component r_3 (2,0)
exp_gasuss1_2_122 = -0.19831503411072457 # matrix component r_2 (3,0)

exp_gasuss1_02 = 0.98366052181486430 # root
exp_gasuss1_2_200 = -0.93754156602517169 # matrix component r_1 (1,0)
exp_gasuss1_2_210 = -0.27632678291615598 # matrix component r_1 (2,0)
exp_gasuss1_2_220 = -9.2911475768203619E-002 # matrix component r_1 (3,0)
exp_gasuss1_2_201 = -0.90268621931130066 # matrix component r_2 (1,0)
exp_gasuss1_2_211 = -0.36816764801502594 # matrix component r_2 (2,0)
exp_gasuss1_2_221 = -0.13114188727343115 # matrix component r_2 (3,0)
exp_gasuss1_2_202 = -0.84968459279566322 # matrix component r_3 (1,0)
exp_gasuss1_2_212 = -0.46907010985682873 # matrix component r_3 (2,0)
exp_gasuss1_2_222 = -0.17345636146458968 # matrix component r_2 (3,0)



def test_calculate_unit_vector():
    """
    A test for the function calculate_unit_vector on the file gauss_method.py, that provided the Right ascension (RA) and 
    Declination (DEC) observation, provided the equatorial heliocentric position of the object.
    Three test are performed each with a different RA and DEC values, to check the return 
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
        
    Due to numpy precise to the 16 digit. 
    """
    result1 = np.array(gm.calculate_unit_vector(test1_ra_dec[0],test1_ra_dec[1]))
    result2 = np.array(gm.calculate_unit_vector(test2_ra_dec[0],test2_ra_dec[1]))
    result3 = np.array(gm.calculate_unit_vector(test3_ra_dec[0],test3_ra_dec[1]))
    result4 = np.array(gm.calculate_unit_vector([test1_ra_dec[0],test2_ra_dec[0],
                                                test3_ra_dec[0]],[test1_ra_dec[1],
                                                test2_ra_dec[1],test3_ra_dec[1]]))
    assert len(result1) == len(result2) == len(result3) == 3
    assert np.array_equal(result1,exp_unit_vec1)
    assert np.array_equal(result2,exp_unit_vec2)
    assert np.array_equal(result3,exp_unit_vec3)
    assert np.array_equal(result4,exp_unit_vec4)

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
    result = gm.calculate_delta_t(test1_times[0],test1_times[1],test1_times[2])
    assert float(result[0][0]) * 1.7202098949957226E-002  == exp_cal_dt1_0 # -6.055580000000191
    assert float(result[1][0]) * 1.7202098949957226E-002 == exp_cal_dt1_1  # 7
    assert round(result[2][0] * 1.7202098949957226E-002,15) == round(exp_cal_dt1_2,15) # 13.055580000000191

def test_calculate_a_and_b():
    """A test for the function calculate_a_and_b on the file gauss_method.py, the tau values of the gauss method it
    generates the necessary a_1, b_1, a_3, and b_3 values, this test is to make the function works correctly 
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
    assert result[0] == exp_cal_ab1_0 
    assert round(result[1][0] * mu** 2,17) == round(exp_cal_ab1_1,17)
    assert round(result[2][0],15) == round(exp_cal_ab1_2,15)
    assert round(result[3][0] *mu **2,17) == round(exp_cal_ab1_3,17)
    
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
    result1 = gm.get_unit_matrix(test1_values_ra_dec[0],test1_values_ra_dec[1])
    assert type(result1) == list 
    assert np.array_equal(np.array(result1).squeeze(),test1_rho)

# Test for function calculate_position are on file test_wis

def test_calculate_values_for_cs():
    """A test for the function calculate_values_for_cs on the file gauss_method.py, This function provided the parameters; 
    the a1, b1, a3, b3, position matrix R, unit matrix of the object, and the "B" matrix 
    provides values necessary to calculate the four coefficient used on the eight degree equation used on the gauss method
    Provided Values:
        test1_a_1 = 0.5361692088746649
        test1_b_1 = 10.85279479419046
        test1_a_3 = 0.4638307911253351
        test1_b_3 = 10.341735205810208
        test1_position_R = np.array([[-0.96969860078090808, -0.94067474928282591 , -0.89451148183020490 ],
                            [-0.22449591050121329, -0.31618137964297799 ,  0.41779882542249258 ],
                            [-9.7312854877537963E-002 ,-0.13705996332878803 ,-0.18111530838131065]])
        test1_rho = np.array([[-0.52583173389865490, -0.58752551582774715, -0.65408632170819792 ],
                    [ 0.84753826708376778 , 0.80401266287853912, 0.74811896536180311 ],
                    [-7.1971337723971657E-002,-9.1528171522763241E-002, -0.11175463041961636]])
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
        than the Fortran Orbfit code.  
            
    """
    result = gm.calculate_values_for_cs(test1_a_1,test1_b_1,test1_a_3,test1_b_3,test1_position_R,test1_rho,test1_B_matrix)
    assert result[0] == exp_cal_valcs1_0
    assert result[1] == exp_cal_valcs1_1
    assert result[2] == exp_cal_valcs1_2
    assert result[3] == exp_cal_valcs1_3

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
    dot product it calculates the four coefficient used on the eight degree equation for the gauss method
    
    Provided Values:
        root 1.0544796939365233
        test1_B_matrix = np.matrix([[115.62989788917604,85.74605554189431,50.08435887453609],
        [-225.65321388125662,-170.38814074716754,-104.27900008385723],
        [ 111.21597772750607,85.55454550094318,54.77127152627216]])
        test1_a_1 = 0.5361692088746649
        test1_b_1 = 10.85279479419046
        test1_a_3 = 0.4638307911253351
        test1_b_3 = 10.341735205810208
        test1_rho = np.array([[-0.52583173389865490, -0.58752551582774715, -0.65408632170819792 ],
                    [ 0.84753826708376778 , 0.80401266287853912, 0.74811896536180311 ],
                    [-7.1971337723971657E-002,-9.1528171522763241E-002, -0.11175463041961636]])
        test1_position_R = np.array([[-0.96969860078090808, -0.94067474928282591 , -0.89451148183020490 ],
                            [-0.22449591050121329, -0.31618137964297799 ,  0.41779882542249258 ],
                            [-9.7312854877537963E-002 ,-0.13705996332878803 ,-0.18111530838131065]])
                            
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
    """
    result1 = gm.calculate_rhos_and_r(test1_1_root,test1_B_matrix,test1_a_1,test1_b_1,test1_a_3,test1_b_3,test1_position_R,test1_rho)
    # print(lenresult1)
    # exit()
    assert len(result1) == 3 and len(result1[1]) == 3 
    assert result1[0] == test1_1_root # root
    assert round(result1[1][0],11)  == round(exp_calc_rhosr1_10,11) # rho 1
    assert round(result1[1][1],11)  == round(exp_calc_rhosr1_11,11)# rho 2
    assert round(result1[1][2],11)  ==  round(exp_calc_rhosr1_12,11)# rho 3
    result1_part2 = result1[2].squeeze()
    assert round(result1_part2[0,0],11)  == round(exp_calc_rhosr1_200,11)# matrix component r_1 (1,0)
    assert round(result1_part2[1,0],12)  == round(exp_calc_rhosr1_210,12)# matrix component r_1 (2,0)
    assert round(result1_part2[2,0],11)  == round(exp_calc_rhosr1_220,11)# matrix component r_1 (3,0)
    assert round(result1_part2[0,1],12)  == round(exp_calc_rhosr1_201,12)# matrix component r_2 (1,0)
    assert round(result1_part2[1,1],12)  == round(exp_calc_rhosr1_211,12)# matrix component r_2 (2,0)
    assert round(result1_part2[2,1],12)  == round(exp_calc_rhosr1_221,12)# matrix component r_2 (3,0)
    assert round(result1_part2[0,2],12)  == round(exp_calc_rhosr1_202,12)# matrix component r_3 (1,0)
    assert round(result1_part2[1,2],12)  == round(exp_calc_rhosr1_212,12)# matrix component r_3 (2,0) # need to fix
    print(result1_part2[1,2])
    assert round(result1_part2[2,2],12) == round(exp_calc_rhosr1_222,12)# matrix component r_2 (3,0)    

def test_gauss_method():
    """A test for the function  on the file gauss_method.py, this function provided three observation times, the rho matrix and position
    of the observer returns the rho and r values using the gauss method (eight degree equation)
    Provided Values:
        test1_times = [58577.489970740738,58583.545550740739,58590.545550740739]
        test1_rho = np.array([[-0.52583173389865490, -0.58752551582774715, -0.65408632170819792 ],
                    [ 0.84753826708376778 , 0.80401266287853912, 0.74811896536180311 ],
                    [-7.1971337723971657E-002,-9.1528171522763241E-002, -0.11175463041961636]])
        test1_position_R = np.array([[-0.96969860078090808, -0.94067474928282591 , -0.89451148183020490 ],
                            [-0.22449591050121329, -0.31618137964297799 ,  -0.41779882542249258 ],
                            [-9.7312854877537963E-002 ,-0.13705996332878803 ,-0.18111530838131065]])
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
    # roots,test_root,r_values,v2,orbit_elements
    assert round(result[0][0],11) == round(exp_gasuss1_00,11)# root
    # print(result[2][0,2,2])
    # exit()
    print(result[0])

    assert round(result[2][0,0,0],11)  == round(exp_gasuss1_2_000,11)# matrix component r_1 (1,0)
    assert round(result[2][0,1,0],12)  == round(exp_gasuss1_2_010,12)# matrix component r_1 (2,0)
    assert round(result[2][0,2,0],11)  == round(exp_gasuss1_2_020,11)# matrix component r_1 (3,0)
    assert round(result[2][0,0,1],12)  == round(exp_gasuss1_2_001,12)# matrix component r_2 (1,0)
    assert round(result[2][0,1,1],11)  == round(exp_gasuss1_2_011,11)# matrix component r_2 (2,0)
    assert round(result[2][0,2,1],12)  == round(exp_gasuss1_2_021,12)# matrix component r_2 (3,0)
    assert round(result[2][0,0,2],12)  == round(exp_gasuss1_2_002,12)# matrix component r_3 (1,0)
    assert round(result[2][0,1,2],12)  == round(exp_gasuss1_2_012,12)# matrix component r_3 (2,0)
    assert round(result[2][0,2,2],12) == round(exp_gasuss1_2_022,12)# matrix component r_2 (3,0)
    
    assert round(result[0][1],11) == round(exp_gasuss1_01,11)# root
    assert round(result[2][1,0,0],11)  == round(exp_gasuss1_2_100,11)# matrix component r_1 (1,0)
    assert round(result[2][1,1,0],12)  == round(exp_gasuss1_2_110,12)# matrix component r_1 (2,0)
    assert round(result[2][1,2,0],11)  == round(exp_gasuss1_2_120,11)# matrix component r_1 (3,0)
    assert round(result[2][1,0,1],12)  == round(exp_gasuss1_2_101,12)# matrix component r_2 (1,0)
    assert round(result[2][1,1,1],12)  == round(exp_gasuss1_2_111,12)# matrix component r_2 (2,0)
    assert round(result[2][1,2,1],12)  == round(exp_gasuss1_2_121,12)# matrix component r_2 (3,0)
    assert round(result[2][1,0,2],12)  == round(exp_gasuss1_2_102,12)# matrix component r_3 (1,0)
    assert round(result[2][1,1,2],12)  == round(exp_gasuss1_2_112,12)# matrix component r_3 (2,0)
    assert round(result[2][1,2,2],12) == round(exp_gasuss1_2_122,12)# matrix component r_2 (3,0)
    
    assert round(result[0][2],11) == round(exp_gasuss1_02,11)# root
    assert round(result[2][2,0,0],11)  == round(exp_gasuss1_2_200,11)# matrix component r_1 (1,0)
    assert round(result[2][2,1,0],12)  == round(exp_gasuss1_2_210,12)# matrix component r_1 (2,0)
    assert round(result[2][2,2,0],11)  == round(exp_gasuss1_2_220,11)# matrix component r_1 (3,0)
    assert round(result[2][2,0,1],12)  == round(exp_gasuss1_2_201,12)# matrix component r_2 (1,0)
    assert round(result[2][2,1,1],12)  == round(exp_gasuss1_2_211,12)# matrix component r_2 (2,0)
    assert round(result[2][2,2,1],12)  == round(exp_gasuss1_2_221,12)# matrix component r_2 (3,0)
    assert round(result[2][2,0,2],11)  == round(exp_gasuss1_2_202,11)# matrix component r_3 (1,0)
    assert round(result[2][2,1,2],12)  == round(exp_gasuss1_2_212,12)# matrix component r_3 (2,0)
    assert round(result[2][2,2,2],12) == round(exp_gasuss1_2_222,12)# matrix component r_2 (3,0)
test_gauss_method()