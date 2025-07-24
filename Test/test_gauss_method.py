import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../Honor_thesis_Code")))
# import gauss_method as gm
from gauss import gauss_method as gm

test1_times = [58577.489970740738,58583.545550740739,58590.545550740739]
def test_calculate_delta_t():
    result = gm.calculate_delta_t(test1_times)
    assert type(result[0]) == float
    assert type(result[1]) == float
    assert type(result[2]) == float
    
    # values not acquired from fortran code
    assert result[0] * 1.7202098949957226E-002  == -0.10416868635938527,15 # -6.055580000000191
    assert result[1] * 1.7202098949957226E-002 == 0.12041469264970059  # 7
    assert round(result[2] * 1.7202098949957226E-002,15) == round(0.22458337900908587,15) # 13.055580000000191
test_calculate_delta_t()

# fitobs result: 
#  Edit: tau13:   0.22458337900908587 
#  Edit: tau1:  -0.10416868635938527     
#  Edit: tau3:   0.12041469264970059  

# python result
# 13.055580000000191      times mue = 0.22458337900908584
# -6.055580000000191      times mue = -0.10416868635938527
# 7.0                     times mue = 0.12041469264970059

# results1  = (-6.055580000000191, 7.0, 13.055580000000191)-0.10416868635938527
0.12041469264970059
0.22458337900908584


test_d_t1 = -6.055580000000191
test_d_t3 = 7
testd_t = 13.055580000000191
mu = 1.7202098949957226E-002
def test_calculate_a_and_b():
    result = gm.calculate_a_and_b(test_d_t1,test_d_t3,testd_t)
    assert result[0] == 0.53616920887466490 
    assert round(result[1] * mu** 2,17) == round(3.2114744736032957E-003,17)
    assert round(result[2],15) == round(0.46383079112533504,15)
    assert round(result[3] *mu **2,17) == round(3.0602457022409226E-003,17)
    
test_calculate_a_and_b()
# fortran results
#  Edit a(1):   0.53616920887466490 
#  Edit b(1):    3.2114744736032957E-003    
#  Edit a(3):   0.46383079112533504 
#  Edit b(3):    3.0602457022409226E-003

# python result:
# a1 = 0.5361692088746649  
# b1 = 10.85279479419046    times mu square = 0.0032114744736032952
# a3 = 0.4638307911253351  
# b3 = 10.341735205810208   times mu square = 0.0030602457022409217


test1_a_1 = 0.5361692088746649
test1_b_1 = 10.85279479419046
test1_a_3 = 0.4638307911253351
test1_b_3 = 10.341735205810208
postion_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
postion_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
postion_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
postion_matrix = np.matrix([[postion_vector1[0],postion_vector2[0],postion_vector3[0]],[postion_vector1[1],postion_vector2[1],postion_vector3[1]],[postion_vector1[2],postion_vector2[2],postion_vector3[2]]])
test1_position_R = postion_matrix.copy()
unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
unit_matrix = np.matrix([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],[unit_vector1[1],unit_vector2[1],unit_vector3[1]],[unit_vector1[2],unit_vector2[2],unit_vector3[2]]])
test1_rho =unit_matrix.copy()
# 115.62989788917604
# 85.74605554189431
# 50.08435887453609
# -225.65321388125662
# -170.38814074716754
# -104.27900008385723
# 111.21597772750607
# 85.55454550094318
# 54.77127152627216
test1_B_matrix = np.matrix([[115.62989788917604,85.74605554189431,50.08435887453609],
 [-225.65321388125662,-170.38814074716754,-104.27900008385723],
 [ 111.21597772750607,85.55454550094318,54.77127152627216]])
def test_calculate_values_for_cs():
    result = gm.calculate_values_for_cs(test1_a_1,test1_b_1,test1_a_3,test1_b_3,test1_position_R,test1_rho,test1_B_matrix)
    assert result[0] == 1.0320244737742499
    assert result[1] == -3527.3938312967603
    assert result[2] == 1.0018109014773862
    assert result[3] == 0.31100143213162895
test_calculate_values_for_cs()

# not proven with fortran code
# results = (1.0320244737742499, -3527.3938312967603, 1.0018109014773862, 0.31100143213162895)
# 1.0320244737742499 -3527.3938312967603 1.0018109014773862 0.31100143213162895




def test_calculate_cs():
    result = gm.calculate_cs(test1_a_1,test1_b_1,test1_a_3,test1_b_3,test1_position_R,test1_rho,test1_B_matrix)
    assert round(-1*result[0],14) == round(-1.0895161396889805,14)
    assert round(-1*result[1],12) == round(2.8036979214270064,12)
    assert round(-1*result[2],12) == round(-2.7106217754654516,12)
    assert round(-1*result[3],14) == 1.0000000000000000

#  Edit c:   -1.0895161396889805     
#  Edit c3:    2.8036979214270064     
#  Edit c6:   -2.7106217754654516     
#  Edit c8:    1.0000000000000000 

# python results
# c = 1.0895161396889808
# c3 =-2.8036979214269007
# c6 = 2.7106217754653157
# c8 = -1


test1_1_sigma = 0.00025237609721519305 # for root 1.0544796939365233
def test_calculate_rhos_and_r():
    result1 = gm.calculate_rhos_and_r(test1_1_sigma,test1_B_matrix,test1_a_1,test1_b_1,test1_a_3,test1_b_3,test1_position_R)
    assert len(result1) == 3 and len(result1[0]) == 1 and type(result1[0][0]) == str
    assert round(result1[1][0],11)  == round(0.13131676742198489,11) # rho 1
    assert round(result1[1][1],11)  == round(0.14179458528999755,11)# rho 2
    assert round(result1[1][2],11)  ==  round(0.15390615730938734,11)# rho 3
    assert round(result1[2][0][0,0],11)  == round(-1.0387491242843767,11)# matrix component r_1 (1,0)
    assert round(result1[2][0][1,0],12)  == round(-0.11319992500134204,12)# matrix component r_1 (2,0)
    assert round(result1[2][0][2,0],11)  == round( -0.10676389829448588,11)# matrix component r_1 (3,0)
    assert round(result1[2][1][0,0],12)  == round(-1.0239826861469132,12)# matrix component r_2 (1,0)
    assert round(result1[2][1][1,0],12)  == round(-0.20217673754220894,12)# matrix component r_2 (2,0)
    assert round(result1[2][1][2,0],12)  == round(-0.15003816245221002,12)# matrix component r_2 (3,0)
    assert round(result1[2][2][0,0],12)  == round(-0.99517939415294532,12)# matrix component r_3 (1,0)
    assert round(result1[2][2][1,0],12)  == round( -0.30265871025338281,12)# matrix component r_3 (2,0)
    assert round(result1[2][2][2,0],12) == round(-0.19831503411072457,12)# matrix component r_2 (3,0)
    
test_calculate_rhos_and_r()
#  fitobs result: 
#  Edit: root:    1.0544796939362728
#  Edit: rho1:   0.13131676742198489
#  Edit: rho2:   0.14179458528999755
#  Edit: rho3:   0.15390615730938734
#  Edit: r1_1:   -1.0387491242843767
#  Edit: r1_2:  -0.11319992500134204
#  Edit: r1_3:  -0.10676389829448588
#  Edit: r2_1:   -1.0239826861469132
#  Edit: r2_2:  -0.20217673754220894
#  Edit: r2_3:  -0.15003816245221002
#  Edit: r3_1:  -0.99517939415294532
#  Edit: r3_2:  -0.30265871025338281
#  Edit: r3_3:  -0.19831503411072457

# python results: 
# root = 1.0544796939365233
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


def test_gauss_method():
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
    #  Edit: root:    1.2457493693119874
#  Edit: rho1:   0.45894799661611524
#  Edit: rho2:   0.49211020370124259
#  Edit: rho3:   0.53267219814382072
#  Edit: r1_1:   -1.2110280216108740
#  Edit: r1_2:   0.16448007923237593
#  Edit: r1_3:  -0.13034395613973659
#  Edit: r2_1:   -1.2298020505564962
#  Edit: r2_2:    7.9481455664558387E-002
#  Edit: r2_3:  -0.18210191046125732
#  Edit: r3_1:   -1.2429250805903169
#  Edit: r3_2:   -1.9296651670140053E-002
#  Edit: r3_3:  -0.24064389301967798
    
    
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
    #  Edit: root:    1.0544796939362728
#  Edit: rho1:   0.13131676742198489
#  Edit: rho2:   0.14179458528999755
#  Edit: rho3:   0.15390615730938734
#  Edit: r1_1:   -1.0387491242843767
#  Edit: r1_2:  -0.11319992500134204
#  Edit: r1_3:  -0.10676389829448588
#  Edit: r2_1:   -1.0239826861469132
#  Edit: r2_2:  -0.20217673754220894
#  Edit: r2_3:  -0.15003816245221002
#  Edit: r3_1:  -0.99517939415294532
#  Edit: r3_2:  -0.30265871025338281
#  Edit: r3_3:  -0.19831503411072457

    
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
    #  Edit: root:   0.98366052181486430
#  Edit: rho1:   -6.1154610272977687E-002
#  Edit: rho2:   -6.4658519414266299E-002
#  Edit: rho3:   -6.8533598008092847E-002
#  Edit: r1_1:  -0.93754156602517169
#  Edit: r1_2:  -0.27632678291615598
#  Edit: r1_3:   -9.2911475768203619E-002
#  Edit: r2_1:  -0.90268621931130066
#  Edit: r2_2:  -0.36816764801502594
#  Edit: r2_3:  -0.13114188727343115
#  Edit: r3_1:  -0.84968459279566322
#  Edit: r3_2:  -0.46907010985682873
#  Edit: r3_3:  -0.17345636146458968
    
    
    # print(result)
#  Edit: root:   0.98366052181486430
#  Edit: rho1:   -6.1154610272977687E-002
#  Edit: rho2:   -6.4658519414266299E-002
#  Edit: rho3:   -6.8533598008092847E-002
#  Edit: r1_1:  -0.93754156602517169
#  Edit: r1_2:  -0.27632678291615598
#  Edit: r1_3:   -9.2911475768203619E-002
#  Edit: r2_1:  -0.90268621931130066
#  Edit: r2_2:  -0.36816764801502594
#  Edit: r2_3:  -0.13114188727343115
#  Edit: r3_1:  -0.84968459279566322
#  Edit: r3_2:  -0.46907010985682873
#  Edit: r3_3:  -0.17345636146458968

#  Edit: root:    1.0544796939362728
#  Edit: rho1:   0.13131676742198489
#  Edit: rho2:   0.14179458528999755
#  Edit: rho3:   0.15390615730938734
#  Edit: r1_1:   -1.0387491242843767
#  Edit: r1_2:  -0.11319992500134204
#  Edit: r1_3:  -0.10676389829448588
#  Edit: r2_1:   -1.0239826861469132
#  Edit: r2_2:  -0.20217673754220894
#  Edit: r2_3:  -0.15003816245221002
#  Edit: r3_1:  -0.99517939415294532
#  Edit: r3_2:  -0.30265871025338281
#  Edit: r3_3:  -0.19831503411072457

#  Edit: root:    1.2457493693119874
#  Edit: rho1:   0.45894799661611524
#  Edit: rho2:   0.49211020370124259
#  Edit: rho3:   0.53267219814382072
#  Edit: r1_1:   -1.2110280216108740
#  Edit: r1_2:   0.16448007923237593
#  Edit: r1_3:  -0.13034395613973659
#  Edit: r2_1:   -1.2298020505564962
#  Edit: r2_2:    7.9481455664558387E-002
#  Edit: r2_3:  -0.18210191046125732
#  Edit: r3_1:   -1.2429250805903169
#  Edit: r3_2:   -1.9296651670140053E-002
#  Edit: r3_3:  -0.24064389301967798

# python results:
# 0.9836605218147869
# -0.06115461027327507
# -0.06465851941455836
# -0.06853359800837579
# -0.9375415660250154
# -0.276326782916408
# -0.09291147576818222
# -0.9026862193111291
# -0.36816764801526075
# -0.13114188727340442
# -0.8496845927954781
# -0.4690701098570404
# -0.17345636146455806

# 1.0544796939365233
# 0.13131676742254747
# 0.1417945852906265
# 0.15390615731009605
# -1.0387491242846727
# -0.11319992500086523
# -0.10676389829452637
# -1.0239826861472827
# -0.20217673754170323
# -0.1500381624522676
# -0.995179394153409
# -0.3026587102528526
# -0.19831503411080378

# 1.2457493693117867
# 0.4589479966158301
# 0.4921102037009649
# 0.5326721981435508
# -1.211028021610724
# 0.16448007923213429
# -0.13034395613971606
# -1.2298020505563332
# 0.07948145566433512
# -0.18210191046123192
# -1.2429250805901404
# -0.019296651670341947
# -0.24064389301964784

