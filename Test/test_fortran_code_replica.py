import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../Honor_thesis_Code")))
import Fortran_Code_Replica as fcr
import numpy as np


unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
test_matrix1 = np.matrix([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],[unit_vector1[1],unit_vector2[1],unit_vector3[1]],[unit_vector1[2],unit_vector2[2],unit_vector3[2]]])
test1_obs_time = [58577.489970740738,58583.545550740739,58590.545550740739]
sqr_mue = 1.7202098949957226E-002

def test_get_tau_values():
    result1 = fcr.get_tau_values(test1_obs_time,sqr_mue)
    assert result1[0] == 0.22458337900908587
    assert result1[1] == -0.10416868635938527
    assert result1[2] == 0.12041469264970059
# test_get_tau_values()
    
#  Edit: tau13:   0.22458337900908587 
#  Edit: tau1:  -0.10416868635938527     
#  Edit: tau3:   0.12041469264970059     

# computer results: 
#  0.22458337900908587
# -0.10416868635938527
# 0.12041469264970059

tau_test1 = [0.22458337900908587,-0.10416868635938527,0.12041469264970059]
def test_make_a_b():
    result1 = fcr.make_a_b(tau_test1[0],tau_test1[1],tau_test1[2])
    assert True
# fortran results
#  Edit a(1):   0.53616920887466490 
#  Edit b(1):    3.2114744736032957E-003    
#  Edit a(3):   0.46383079112533504 
#  Edit b(3):    3.0602457022409226E-003

# python results
# a1  = 0.5361692088746649
# b_1 = 0.0032114744736032957 
# a_3 = 0.46383079112533504
# b_3 = 0.0030602457022409226


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




# -181.39997132286484   346.8208661801217   -167.22601164748204
# -49.141453222161964   99.18166675238325   -49.58309521064877
# 732.7444648345615   -1365.948759032001   637.8813360694052
test1_inv_matr = [[-181.39997132286484,-49.141453222161964,732.7444648345615],
                [346.8208661801217,99.18166675238325,-1365.948759032001],
                [-167.22601164748204,-49.58309521064877,637.8813360694052]]

test1_dot_b =[[-0.005851587221372763],
[-0.001999529945903465],
[-0.0008667750934766205]]
test1_dot_a = [[0.005850249367257254],
[0.002025625186994423],
[0.0008769501443510741]]
test1_obs_postions = matrix1_test3.copy()
test1_unit_vector = test_matrix1.copy()
def test_get_d_and_pos_variables():
    result1 = fcr.get_d_and_pos_variables(np.matrix(test1_inv_matr),np.matrix(test1_dot_b),np.matrix(test1_dot_a),np.matrix(test1_obs_postions),np.matrix(test1_unit_vector))
    assert result1[0] == 1.0320244737743007
    assert result1[1] == -1.0437988981068052
    assert result1[2] == 1.0036250823189330 
    assert result1[3] == 0.31100143213162895
    
# results from fortran
#  Edit: a2star:    1.0320244737743007    
#  Edit: b2star:   -1.0437988981068052     
#  Edit: r22:    1.0036250823189330     
#  Edit: s2r2:   0.31100143213162895 

# python results:
# a2star = 1.0320244737743007
# b2star = -1.0437988981068052
# r22  = 1.003625082318933
# s2r2 = 0.31100143213162895


test1_d_1 =1.0320244737743007
test1_d_2 = -1.0437988981068052
test1_r2pos = 1.003625082318933 
test1_dotr2 = 0.31100143213162895

def test_generate_Cs():
    result1 = fcr.generate_Cs(test1_d_1,test1_d_2,test1_r2pos,test1_dotr2)
    assert len(result1) == 4
    assert result1[0] == 1.0000000000000000 
    assert result1[1] == -2.7106217754654516
    assert result1[2] == 2.8036979214270064
    assert result1[3] == -1.0895161396889805
# test_generate_Cs()

# fortran results:
#  Edit c:   -1.0895161396889805     
#  Edit c3:    2.8036979214270064     
#  Edit c6:   -2.7106217754654516     
#  Edit c8:    1.0000000000000000 

# python results
# -1.0895161396889805
# 2.8036979214270064
# -2.7106217754654516
# 1

test1_c = -1.0895161396889805
test1_c3 =2.8036979214270064
test1_c6 =-2.7106217754654516
test1_c8 =1


# test root 1.054479693936294
test1_a1 = 0.5361692088746649
test1_b1 = 0.0032114744736032957
test1_a3 = 0.46383079112533504
test1_b3 = 0.0030602457022409226
test1_r2m3 = 0.8528749073206177
test1_pos_matrix = matrix1_test3.copy()
test1_inv_matrix = test1_inv_matr.copy()
def test_calculate_rhos():
    result1 = fcr.calculate_rhos(test1_a1,test1_b1,test1_a3,test1_b3,test1_r2m3,test1_pos_matrix, test1_inv_matrix)
    # (float(p_1),float(p_3),c_1,c_3,gcap,crhom,r2m3)
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
# fortran results
#  Edit root:    1.0544796939362728     
#  Edit c_1:   0.53890819486870201     
#  Edit c_3:   0.46644079789501225     
#  Edit r2m3:   0.85287490732066884     
#  Edit gcap:    8.5957745815007369E-004   3.2027626969707845E-004   1.3769941683434983E-004
#  Edit crhom:   -7.0767682087375050E-002  0.14179458528999755       -7.1788110816345896E-002
#  Edit rho1:   0.13131676742198489     
#  Edit rho3:   0.15390615730938734 

# python results
# 1.054479693936294
# 0.5389081948687019
# 0.4664407978950121
# 0.8528749073206177
# 0.0008595774581503512     0.0003202762696971895    0.00013769941683439146
# -0.07076768208740036      0.14179458529004793     -0.07178811081637124
# 0.13131676742203188
# 0.1539061573094417


test_roots1 = 1.245749369311985
test_a1 = 0.5361692088746649
test_b1 = 0.0032114744736032957
test_a3 = 0.46383079112533504
test_b3 = 0.0030602457022409226
test_pos_matrix = matrix1_test3.copy()
test_inv_matrix = test1_inv_matr.copy()
test1_unit_vector = test_matrix1.copy()
def test_find_rho_and_r():
    result1 = fcr.find_rho_and_r(test_roots1,test_a1,test_b1,test_a3,test_b3,np.matrix(test_pos_matrix),test_inv_matrix,test1_unit_vector)    #return(rhos,r)
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
    # assert type

# test_find_rho_and_r()

# fortran results
#  Edit: root:    1.2457493693119874     
#  Edit: rho1:   0.45894799661611524     
#  Edit: rho3:   0.53267219814382072     
#  Edit: r1_1:   -1.2110280216108740     
#  Edit: r1_2:   0.16448007923237593     
#  Edit: r1_3:  -0.13034395613973659     
#  Edit: r3_1:   -1.2429250805903169     
#  Edit: r3_2:   -1.9296651670140053E-002
#  Edit: r3_3:  -0.24064389301967798   
# python results
# rho1
# 0.45894799661611524
# rho3
# 0.5326721981438207
# r_1
# -1.211028021610874
# 0.16448007923237593
# -0.1303439561397366
# r_3
# -1.242925080590317
# -0.019296651670140053
# -0.24064389301967798

def test_gauss_method_8th():
    results1 = fcr.gauss_method_8th(np.matrix(test_pos_matrix),test1_obs_time,np.matrix(test1_unit_vector))
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
    
# Edit: root:    1.2457493693119874 
# Edit: root:    1.0544796939362728  
# Edit: root:   0.98366052181486430  
# rho values
# Edit: rho1:   0.45894799661611524     
#  Edit: rho3:   0.53267219814382072 
# Edit: rho1:   0.13131676742198489     
#  Edit: rho3:   0.15390615730938734   
#  Edit: rho1:   -6.1154610272977687E-002
#  Edit: rho3:   -6.8533598008092847E-002
# r values
# r1
# Edit: r1_1:   -1.2110280216108740     
#  Edit: r1_2:   0.16448007923237593     
#  Edit: r1_3:  -0.13034395613973659  
# Edit: r1_1:   -1.0387491242843767     
#  Edit: r1_2:  -0.11319992500134204     
#  Edit: r1_3:  -0.10676389829448588   
#  Edit: r1_1:  -0.93754156602517169     
#  Edit: r1_2:  -0.27632678291615598     
#  Edit: r1_3:   -9.2911475768203619E-002
#  r3
#  Edit: r3_1:   -1.2429250805903169     
#  Edit: r3_2:   -1.9296651670140053E-002
#  Edit: r3_3:  -0.24064389301967798 
#  Edit: r3_1:  -0.99517939415294532     
#  Edit: r3_2:  -0.30265871025338281     
#  Edit: r3_3:  -0.19831503411072457   
#  Edit: r3_1:  -0.84968459279566322     
#  Edit: r3_2:  -0.46907010985682873     
#  Edit: r3_3:  -0.17345636146458968     
    
# Python result
# 1.245749369311985
# 1.054479693936294
# 0.9836605218148519
# 0.45894799661611524
# 0.5326721981438207
# 0.13131676742203188
# 0.1539061573094417
# -0.061154610273003375
# -0.06853359800812431
# -1.211028021610874
# 0.16448007923237593
# -0.1303439561397366
# -1.0387491242844016
# -0.11319992500130222
# -0.10676389829448926
# -0.9375415660251583
# -0.27632678291617774
# -0.09291147576820177
# -1.242925080590317
# -0.019296651670140053
# -0.24064389301967798
# -0.995179394152981
# -0.3026587102533421
# -0.19831503411073065
# -0.8496845927956426
# -0.46907010985685227
# -0.17345636146458615
