'''This file contains the functions utlized to perform the eight degree equation of 
    the gauss method, and also includes the calculation fo the orbital elments 
    based on the calculated roots (file is still in progress).
'''

import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from gauss import math_fucntions as mf
from astropy.time import Time
module_path = "/home/mila/mpc-software/wis/wis" # contains path for wis module
if module_path not in sys.path:
    sys.path.append(module_path)
np.set_printoptions(threshold=np.inf,precision=20)

def mjd_to_jd(date):
    """convers a modified julian date to, standart julian date
    Args:
        date (float): The date to be converted 

    Returns:
        date: the date in Julian date
    """
    if type(date) != list and type(date) != np.ndarray:
        date = [date]
    date1 = np.array(date)
    date1 = date1 + 2400000.5
    return np.array(date1)


def calculate_position_obs(mjd,obscode): # uses Wis.py
    """Uses Wis to calculate the position of the obsevatory in a Heliocentric equatorial frame

    Args:
        mjd (list/int): list of dates or a single date 
        obscode (string): Observatory code

    Returns:
        np.matrix: a 3x3 matrix per each of the given times 
        (such as provided a list of 4 times, it returns a list with 4, 3x3 matrices corresposing to each time)
    """
    # import time
    # start_time = time.time()
    jd_times = mjd_to_jd(mjd)
    times = Time(jd_times, format='jd', scale='tdb',)
    import wis
    W = wis.wis(str(obscode), times)
    observer_posns = W.hXYZ
    # end_time = time.time()
    # execution_time = end_time - start_time
    return np.matrix(observer_posns)
# x = calculate_position_obs([51975.54306,54275.54306,54930.54306,54175.54306,54276.54306,54475.54306],500) # or a single time
# x = calculate_position_obs([51975.54306,54275.54306],500) # or a single time
    


def calculate_unit_vector(ra,dec):
    """Provided the declination and right ascession calculates the unti vector

    Args:
        ra (float): Right Ascension
        dec (float): Declination

    Returns:
        array: array of lenght three with the unit vector of provided RA and DEC
    """
    x_component = np.cos(ra) * np.cos(dec) #x hat
    y_component = np.sin(ra) * np.cos(dec) #y hat
    z_component = np.sin(dec) #z hat
    xyz_unit_vector = np.array([x_component,y_component,z_component]).T
    return xyz_unit_vector

def get_unit_matrix(values):
    """provided a list of RA and dec Values it generates a matrix 
    it uses function calculate_unit_vector()

    Args:
        values (List):list where each elements is a sublist with ra and dec values, in the follwoing format:
        [[Ra1,Dec1],[Ra2,Dec2]...] 

    Returns:
        np,matrix: a matrix containng the unit vectors of the provided coordinates
    """
    values_t = np.array(values).T
    if len(values_t[0]) % 3 != 0 or len(values_t[0]) != len(values_t[1]):
        print("The provided RA and DEC values are not compatible")
        print(values)
        exit()
    ra = values_t[0]
    dec = values_t[1]
    unit_matrices = calculate_unit_vector(ra,dec)
    if len(unit_matrices) != len(values):
        print("The provided RA and DEC values are not compatible 1")
        print(values)
        exit()
    else:
        try:
            unit_matrices = np.transpose(unit_matrices, axes=(0, 2, 1))
        except ValueError:
            unit_matrices = unit_matrices.T
    return unit_matrices
    
test21_ra_dec = [1.1536801615828238,0.20914571314856717]
test22_ra_dec = [1.2043302709741688,0.21933455747076516]
test23_ra_dec = [1.2623223260301843,0.22963781782170506]
test2_pos = [58577.489970740738,58583.545550740739,58590.545550740739]
test1_ra_dec = [2.1260971174823160,-7.2033616739254860E-002]
test2_ra_dec =[2.2018577151544165,-9.1656450482163324E-002]
test3_ra_dec = [2.2892342062190609,-0.11198856664053504]
test1_values_ra_dec = [test1_ra_dec,test2_ra_dec,test3_ra_dec]
test2_values_ra_dec = [test21_ra_dec,test22_ra_dec,test23_ra_dec]

mue = 1.7202098949957226E-002 **2 # Only global variable

def calculate_delta_t(times):
    """Calculates the detla t's values of the gauss method

    Args:
        times (list): a list of 3 times with the following format [times1, times2,times3]

    Returns:
        tuple: a list of the tau values.
    """
    times = np.array(times)
    if type(times[0]) == list or type(times[0]) == np.ndarray:
        times0 = times[:,0]
        times1 = times[:,1]
        times2 = times[:,2]
    else:
        times0 = times[0]
        times1 = times[1]
        times2 = times[2]
    d_t1 =times0 - times1 
    d_t3 = times2 - times1
    d_t = times2 - times0 
    return (d_t1,d_t3,d_t)

times1 = [58577.489970740738,58583.545550740739,58590.545550740739]
times2 = [58577.489970740738,58601.545550740739,58626.545550740739]
# calculate_delta_t([times1,times2])

def calculate_a_and_b(d_t1,d_t3,d_t):
    """calculates a and b values, of the gauss method using it's appropiate delta t values

    Args:
        d_t1 (float): delta one value
        d_t3 (float): delta three value
        d_t (float): tau or delta t value

    Returns:
        tupple: tupple with the a's and b's values.
    """
    d_t1 = np.array(d_t1)
    d_t3 = np.array(d_t3)
    d_t = np.array(d_t)
    a_1 = d_t3 / d_t
    b_1 = (d_t3 * (d_t**2 - d_t3**2)) / (6 * d_t)
    a_3 = - d_t1 / d_t
    b_3 = - (d_t1 * (d_t**2 - d_t1 **2)) / (6 * d_t)
    return (a_1,b_1,a_3,b_3)

d_t = 13.055580000000191,49.05558000000019
d_t1 = -6.055580000000191, -24.05558000000019
d_t3 = 7.0,25.0
# calculate_a_and_b(d_t1,d_t3,d_t)

def calculate_values_for_cs(a_1,b_1,a_3,b_3,position_R,rho,B_matrix):
    """caculate values used to genrate the c values for the gauss method (c1,c3,c6,c8) 
    Args:
        a_1 (float): a one value
        b_1 (float): b one value
        a_3 (float): a three value
        b_3 (float): b three value
        position_R (np.matrix): The position of the observer 
        rho (np.matrix): The unit vector of the matrix
        B_matrix (np.matrix): matrix  used for the calculation

    Returns:
        tupple: d_1 and d_2 values used to generate c values, and also additional retunrs for cheking oon test 
    """
    
    if type(a_1) == np.ndarray or type(a_1) == list and len(a_1) == len(rho):
        a_1 = np.array(a_1)
        b_1 = np.array(b_1)
        a_3 = np.array(a_3)
        b_3 = np.array(b_3)
        position_R = np.array(position_R)
        rho = np.array(rho)
        B_matrix = np.array(B_matrix)
        rho = np.array(rho,dtype = np.matrix)
        rh01 = np.array(rho[:,0,1])
        rh11 = np.array(rho[:,1,1])
        rh21 = np.array(rho[:,2,1])
        B_matrix = np.array(B_matrix,dtype = np.matrix)
        bmatrix10 = np.array(B_matrix[:,1,0])
        bmatrix11 = np.array(B_matrix[:,1,1])
        bmatrix12 = np.array(B_matrix[:,1,2])
        position_R = np.array(position_R,dtype = np.matrix)
        pos_r01 = np.array(position_R[:,0,1])
        pos_r11 = np.array(position_R[:,1,1])
        pos_r21 = np.array(position_R[:,2,1])
        temp_matrix1 = np.array([[rh01],[rh11],[rh21]],dtype = np.matrix).T # vector components of p_2
        temp_matrix2 = np.array([[pos_r01,pos_r11,pos_r21]],dtype = np.matrix).T # vecotr componets of R_2
    else: 
        a_1 = np.array(a_1)
        b_1 = np.array(b_1)
        a_3 = np.array(a_3)
        b_3 = np.array(b_3)
        position_R = np.array(position_R)
        rho = np.array(rho)
        B_matrix = np.array(B_matrix)
        rh01 = rho[0,1]
        rh11 = rho[1,1]
        rh21 =rho[2,1]
        pos_r01 = position_R[0,1]
        pos_r11 = position_R[1,1]
        pos_r21 = position_R[2,1]
        bmatrix10 = B_matrix[1,0]
        bmatrix11 = B_matrix[1,1]
        bmatrix12 = B_matrix[1,2]
        temp_matrix1 = np.matrix([[rh01,rh11,rh21]]) # vector components of p_2
        temp_matrix2 = np.matrix([[pos_r01],[pos_r11],[pos_r21]])# vecotr componets of R_2
        
    d_1 = bmatrix10 * a_1 - bmatrix11 + bmatrix12 * a_3
    d_2 = bmatrix10 * b_1 + bmatrix12 * b_3
    part1_mag = pos_r01** 2 + pos_r11** 2 + pos_r21** 2
    part1_mag = part1_mag.astype(float)
    magnitude_R_2 = np.sqrt(part1_mag)
    pos_2_dot_p2 = temp_matrix1 @ temp_matrix2
    pos_2_dot_p2 = pos_2_dot_p2.squeeze()
    return (d_1,d_2,magnitude_R_2,pos_2_dot_p2)

onea1 =  0.5361692088746649
oneb1 = 10.85279479419046
onea3 =   0.4638307911253351
oneb3 = 10.341735205810208
b1_matrix = np.matrix([[ 115.62989798187795,   85.74605563851054,   50.08435897098997],
                       [-225.65321405248952, -170.38814092604792, -104.27900026288157],
                       [ 111.2159778069043,    85.5545455841151,    54.77127160975296]])
r1 = np.matrix([[-0.9696986008561639,  -0.9406747493915482,  -0.8945114819710074 ],
                [-0.2244959102068138,  -0.31618137935660046, -0.4177988251583791 ],
                [-0.09731285474991118, -0.13705996320464234, -0.18111530826682148]])
rho1 = np.matrix([[-0.5258317338986549,  -0.5875255158277471,  -0.6540863217081979 ],
        [ 0.8475382670837678,   0.8040126628785391,   0.7481189653618031 ],
        [-0.07197133772397166, -0.09152817152276324, -0.11175463041961636]])
# result = calculate_values_for_cs(onea1,oneb1,onea3,oneb3,r1,rho1,b1_matrix)
# print(result)

twoa1 = 0.5096260201184025
twob1 = 151.3122062376672
twoa3 = 0.49037397988159753
twob3 = 149.3825437623352
b2_matrix = np.matrix([[ 24.029891036285346,  20.895058685538306,  13.983555837347513],
             [-48.134551467186526, -40.81300622688785,  -26.069774688569154],
             [ 24.13169102456899,   19.60430664019454,   11.474054536074508]])
r2 = np.matrix([[-0.9696986008561639,  -0.796540944160929,   -0.4773994619640845 ],
      [-0.2244959102068138,  -0.564914570091051,   -0.8191737144498188 ],
      [-0.09731285474991118, -0.24488569907572702, -0.3551090002088778 ]])
rho2 = np.matrix([[0.646743628403743,   0.5248759921239102,  0.3627509124021436 ],
                  [0.6000699731114835,  0.7060661103062899,  0.8015212662825536 ],
                  [0.47078520207111946, 0.4753691626187874,  0.4753687360862347 ]])
# result = calculate_values_for_cs(twoa1,twob1,twoa3,twob3,r2,rho2,b2_matrix)
# print(result)
# vv = calculate_values_for_cs([onea1,twoa1],[oneb1,twob1],[onea3,twoa3],[oneb3,twob3],[r1,r2],[rho1,rho2],[b1_matrix,b2_matrix])
# print(result)

def calculate_cs(d_1,d_2,magnitude_R_2,pos_2_dot_p2):
    """Generates non-zero c values for the Gauss method c,c3,c6,c8 
    uses function calculate_values_for_cs

    Args:
        a_1 (float): a one value
        b_1 (float): b one value
        a_3 (float): a three value
        b_3 (float): b three value
        position_R (np.matrix): The position of the observer 
        rho (np.matrix): The unit vector of the matrix
        B_matrix (np.matrix): matrix  used for the calculation
        
    Returns:
        tupple: with generated c values 
    """
    d_1 = np.array(d_1)
    d_2 = np.array(d_2)
    magnitude_R_2 = np.array(magnitude_R_2)
    pos_2_dot_p2 = np.array(pos_2_dot_p2)
    c = d_2 ** 2 * mue**2
    c_3 =mue * 2 * (d_1 * d_2 + d_2 * pos_2_dot_p2)
    c_6 = d_1**2 +magnitude_R_2**2 + 2 * d_1 *pos_2_dot_p2
    c_8 = -1    # it can be -1 and, or multiply all the other coefficients by -1 instead.
    c = c.squeeze()
    c_3 = c_3.squeeze()
    c_6 = c_6.squeeze()
    return(c,c_3,c_6,c_8)

onedd1 = 1.0320244778078234
onedd2 = -3527.393835006538
onemag = 1.0018109014721053
onepos2 = 0.31100143241439426
# x = calculate_cs(onedd1,onedd2,onemag,onepos2)
twodd1 = 3.498447163830976
twodd2 = -11177.714437049934
twomag = 1.0067645965357288
twopos2 = -0.9333633612002668
# x = calculate_cs(twodd1,twodd2,twomag,twopos2)
# x = calculate_cs([onedd1,twodd1],[onedd2,twodd2],[onemag,twomag],[onepos2,twopos2])

def calculate_rhos_and_r(root,B_matrix,a_1,b_1,a_3,b_3,R,p):
    """provided the appropite values calculates rhos and r values.

    Args:
        sigma (float): calculated using the current root in question (r_2)
        B_matrix (np.matrix): matrix  used for the calculation
        a_1 (float): a one value
        b_1 (float): b one value
        a_3 (float): a three value
        b_3 (float): b three value
        R (np.matrix): The position of the observer 
        P (np.matrix): The unit vector of the matrix
    Returns:
        list: list with 3 of lenght the firt elemnt is a place holder, 
        second element contains all the calculated rho values, 3 elemnts 
        contains calculated r values, where each r value is a 3x3 matrix
    """
    if (type(root) == list or type(root) == np.ndarray) and len(a_3) == len(b_3):
        sigma = np.array(mue/(np.array(root))**3,dtype=float).squeeze()
        a_1 = np.array(a_1)
        b_1 = np.array(b_1)
        a_3 = np.array(a_3)
        b_3 = np.array(b_3)
        B_matrix = np.array(B_matrix,dtype=np.matrix).squeeze()
        b_matrix00 = B_matrix[:,0,0]
        b_matrix10 = B_matrix[:,1,0]
        b_matrix20 = B_matrix[:,2,0]
        b_matrix01 = B_matrix[:,0,1]
        b_matrix11 = B_matrix[:,1,1]
        b_matrix21 = B_matrix[:,2,1]
        b_matrix02 = B_matrix[:,0,2]
        b_matrix12 = B_matrix[:,1,2]
        b_matrix22 = B_matrix[:,2,2] 
        R = np.array(R,dtype=np.matrix).squeeze()
        r0 = R[:,:,0]
        r1 =  R[:,:,1]
        r2 =  R[:,:,2]
        p = np.array(p,dtype=np.matrix).squeeze()
        p0 =  p[:,:,0]
        p1 = p[:,:,1]
        p2 = p[:,:,2]
        r = np.zeros((len(sigma),3,3))
        num = 1
    else:
        num = 0
        sigma = mue/(np.array(root))**3
        r = np.zeros((1,3,3))
        b_matrix00 = B_matrix[0,0]
        b_matrix10 = B_matrix[1,0]
        b_matrix20 = B_matrix[2,0]
        b_matrix01 = B_matrix[0,1]
        b_matrix11 = B_matrix[1,1]
        b_matrix21 = B_matrix[2,1]
        b_matrix02 = B_matrix[0,2]
        b_matrix12 = B_matrix[1,2]
        b_matrix22 = B_matrix[2,2] 
        r0 = np.array(R[:,0]).T
        r1 =  np.array(R[:,1]).T
        r2 =  np.array(R[:,2]).T
        p0 =  np.array(p[:,0])
        p1 = np.array(p[:,1])
        p2 = np.array(p[:,2])  
    mag_p1 = (b_matrix00 * a_1 - b_matrix01 + b_matrix02 * a_3 + (b_matrix00 * b_1 + b_matrix02 * b_3) * sigma)/((-1)*(a_1+b_1*sigma))
    mag_p2 = (b_matrix10 * a_1 - b_matrix11 + b_matrix12 * a_3 + (b_matrix10 * b_1 + b_matrix12 * b_3) * sigma)
    mag_p3 = (b_matrix20 * a_1 - b_matrix21 + b_matrix22 * a_3 + (b_matrix20 * b_1 + b_matrix22 * b_3) * sigma)/ ((-1)* (a_3+b_3*sigma))
    r[:,:,0] = np.array(r0) + np.array(np.matrix(mag_p1).T)*np.array(p0)
    r[:,:,1] = np.array(r1) + np.array(np.matrix(mag_p2).T)*np.array(p1)
    r[:,:,2] = np.array(r2) + np.array(np.matrix(mag_p3).T)*np.array(p2)
    magnitudes = np.array([mag_p1,mag_p2,mag_p3])
    if num == 0:
        r = r.squeeze()
    return (root, magnitudes,r)

s1 = 0.0003109046030227718
b1_ma = np.matrix([[ 115.62989798187795,   85.74605563851054 ,  50.08435897098997],
         [-225.65321405248952, -170.38814092604792, -104.27900026288157],
         [ 111.2159778069043  ,  85.5545455841151,    54.77127160975296]])
oneaa1 = 0.5361692088746649
onebb1 = 10.85279479419046
oneaa3 = 0.4638307911253351
onebb3 = 10.341735205810208
one_r = np.matrix([[-0.9696986008561639,  -0.9406747493915482,  -0.8945114819710074 ],
         [-0.2244959102068138 , -0.31618137935660046, -0.4177988251583791 ],
         [-0.09731285474991118, -0.13705996320464234 ,-0.18111530826682148]])
one_p = np.array([[-0.5258317338986549,  -0.5875255158277471 , -0.6540863217081979 ],
        [ 0.8475382670837678  , 0.8040126628785391,   0.7481189653618031 ],
        [-0.07197133772397166 ,-0.09152817152276324, -0.11175463041961636]])
# calculate_rhos_and_r(s1,b1_ma,oneaa1,onebb1,oneaa3,onebb3,one_r,one_p)

s2 = 0.0003107646078995292
b2_ma = np.matrix([[ 24.029891036285346 , 20.895058685538306  ,13.983555837347513],
                   [-48.134551467186526 ,-40.81300622688785  ,-26.069774688569154],
                   [ 24.13169102456899 ,  19.60430664019454 ,  11.474054536074508]])
twoaa1 = 0.5096260201184025
twobb1 = 151.3122062376672
twoaa3 = 0.49037397988159753
twobb3 = 149.3825437623352
two_r = np.matrix([[-0.9696986008561639 , -0.796540944160929  , -0.4773994619640845 ],
                   [-0.2244959102068138  ,-0.564914570091051   ,-0.8191737144498188 ],
                   [-0.09731285474991118 ,-0.24488569907572702, -0.3551090002088778 ]])
two_p = np.array([[0.646743628403743   ,0.5248759921239102,  0.3627509124021436 ],
                  [0.6000699731114835  ,0.7060661103062899 , 0.8015212662825536 ],
                  [0.47078520207111946, 0.4753691626187874 , 0.4753687360862347 ]])
# calculate_rhos_and_r(s2,b2_ma,twoaa1,twobb1,twoaa3,twobb3,two_r,two_p)

# calculate_rhos_and_r([s1,s2],[b1_ma,b2_ma],[oneaa1,twoaa1],[onebb1,twobb1],[oneaa3,twoaa3],[onebb3,twobb3],[one_r,two_r],[one_p,two_p])


def Check_root(info_root,spurious_distance= 0.001):
    """Checks for SPurious roots with Spurious distance being the 
    earth Hill's sphere

    Args:
        info_root (list): Contains the information of the root in a list as:
        [root,[rho1,rho2,rho3],[matrix r1, matrixr2, matrix r3]]
        spurious_distance (float, optional): the Eartch Hill's sphere in AU. Defaults to 0.001.

    Returns:
        int: if it return one the root is spurious (and should be reomoved), otherwise the root is not spurious. 
    """
    if type(info_root[0]) == list or type(info_root[0]) == np.ndarray:
        # info_root = np.array(info_root)
        # root = np.array([info[0] for info in info_root])    # had to use for loop due to structure of info_root, cannot become a np.array
        # print(root)
        # print(type(info_root))
        # print(root)
        rho_2 = np.array([info[1][1] for info in info_root])
        root = []
        rho_2 = []
        for info in info_root:
            root.append(info[0])
            rho_2.append(info[1][1])
        root = np.array(root)
        rho_2 = np.array(rho_2)
    else:
        root = info_root[0]
        rho_2 = info_root[1][1]
    bool_value = rho_2 < spurious_distance
    # false keep the root true remove the root
    return bool_value

info1 = [1.2457493837871696, [np.float64(0.45894801772552885), np.float64(0.4921102259878265), np.float64(0.5326722220103743)], np.array([[-1.2110280327861294  , -1.2298020637591551  ,
        -1.2429250963419056  ],
       [ 0.16448009741781122 ,  0.07948147386963161 ,
        -0.019296633551005138],
       [-0.13034395753138256 , -0.1821019123769619  ,
        -0.2406438955723867  ]])]
# Check_root(info1)

info2 = [2.3364268105073442, [np.float64(3.2326778980657243), np.float64(3.239113228824863), np.float64(3.1490124057944495)], np.array([[1.1210152323994478, 0.9035918254202033, 0.6649076614035214],
       [1.7153370291635723, 1.7221135082269676, 1.7048266965820191],
       [1.4245840627218025, 1.294888844138185 , 1.141833047053503 ]])]
# Check_root(info2)

# Check_root([info1,info2])



def provide_position_matrix():
    postion_vector1 = [-0.96969860078090808,-0.22449591050121329,-9.7312854877537963E-002]
    postion_vector2 = [-0.94067474928282591,-0.31618137964297799,-0.13705996332878803]
    postion_vector3 = [-0.89451148183020490,-0.41779882542249258,-0.18111530838131065]
    postion_matrix = np.matrix([[postion_vector1[0],postion_vector2[0],postion_vector3[0]],[postion_vector1[1],postion_vector2[1],postion_vector3[1]],[postion_vector1[2],postion_vector2[2],postion_vector3[2]]])
    return postion_matrix

def provide_unit_matrix():
    unit_vector1 = [-0.52583173389865490,0.84753826708376778,-7.1971337723971657E-002]
    unit_vector2 = [-0.58752551582774715,0.80401266287853912,-9.1528171522763241E-002]
    unit_vector3 = [-0.65408632170819792,0.74811896536180311,-0.11175463041961636]
    unit_matrix = np.array([[unit_vector1[0],unit_vector2[0],unit_vector3[0]],[unit_vector1[1],unit_vector2[1],unit_vector3[1]],[unit_vector1[2],unit_vector2[2],unit_vector3[2]]])
    return unit_matrix

def calculated_v2(position_matrix,t1,t3):
    # uses gibbs method to calculate the v2 values
    if type(t1) == list or type(t1) == np.ndarray and type(t3) == list or type(t3) == np.ndarray and len(t1) == len(t3):
        t1 = np.array(t1)
        t3 = np.array(t3)
        r = np.array(position_matrix.copy()).squeeze()
        one = False
         
    else:
        r = np.array([position_matrix.copy()])          
        one = True
    t1 = np.sqrt(mue) * t1
    t3 = np.sqrt(mue) * t3
    t13 = t3 - t1
    normr = []
    for matrix in r:
        matrix = np.array(matrix,dtype=float)
        normr.append( np.array(1/(np.array([np.linalg.norm(matrix,axis=0)]))**3).squeeze())
    normr = np.array(normr)
    if not one:
        normr = normr.squeeze()
    r1m3 = normr[:,0]
    r2m3 = normr[:,1]
    r3m3 = normr[:,2]
    d1= t3 * (r1m3 / 12 - 1/ (t1 * t13)) 
    d2= (t1+t3) * (r2m3 / 12 - 1/ (t1 * t3)) 
    d3= -t1 * (r3m3 / 12 + 1/ (t3 * t13)) 
    v2 =[]
    for i in range(3) :
        v2.append(np.sqrt(mue)*(-d1*r[:,i,0]+d2*r[:,i,1]+d3*r[:,i,2]))
            #   v2(i)=gk*(-d1*r(i,1)+d2*r(i,2)+d3*r(i,3)) 
    if not one:
        v2 = np.array(v2).T
    v2 = np.array(v2).squeeze()
    return v2

d1t1 = -6.055580000000191
d1t3 = 7.0
pos1 = [[-1.2110280327861294  , -1.2298020637591551 ,  -1.2429250963419056  ],
        [ 0.16448009741781122  , 0.07948147386963161 , -0.019296633551005138],
        [-0.13034395753138256  ,-0.1821019123769619   ,-0.2406438955723867  ]]
# calculated_v2(pos1,d1t1,d1t3)

d2t1 = -6.055580000000191
d2t3 = 7.0
pos2 = [[-1.0387491066359327,  -1.0239826649099535,  -0.9951793684967545 ],   
        [-0.1131999532740922 , -0.2021767664668157  ,-0.30265873949488553],
        [-0.10676389574099107 ,-0.15003815900270875 ,-0.198315029588662  ]]
# calculated_v2(pos2,d2t1,d2t3)

# calculated_v2([pos1,pos2],[d2t1,d1t1],[d2t3,d1t3])

# calculation the orbital elments. 
eclip = np.radians(23.439281)# used for coordinate rotation

def rotate_from_equeritorial_to_ecliptic(vect):
    if type(vect[0]) != float and type(vect[0]) != int and (type(vect[0]) == list or type(vect[0]) == np.ndarray):
        copies = len(vect)
    else:
        copies = 1
    rot_matrix = np.array([[1,0,0],[0,np.cos(-eclip),-1*np.sin(-eclip)],[0,np.sin(-eclip),np.cos(-eclip)]])
    rotation = np.array([np.copy(rot_matrix) for _ in range(copies)])
    rotation = rotation.squeeze()
    vect = np.array(vect)
    vect = np.reshape(vect, (copies, 3,1))
    res = rotation @ vect
    res = np.reshape(res, (copies, 3))
    res = res.squeeze()
    return res



ve1 = [-1.2298020637591551  , 0.07948147386963161 ,-0.1821019123769619 ]
ve2 = [-0.0025254932789341334, -0.014086836069343118 , -0.00846957995751123  ]
ve3 = [ 0.0032252372826048928 ,-0.014561432694412485  ,-0.007042477567136881 ]
# rotate_from_equeritorial_to_ecliptic(ve1)
# rotate_from_equeritorial_to_ecliptic([ve1,ve2])

def semi_major_axis(r2,v2):
    # only accurate for elliptic orbits, not hyperbolic orbits
    r2 = np.array(r2)
    v2 = np.array(v2)
    a = ((2/r2)-(v2**2/mue))**(-1) # units AU
    return a

or2 = 1.2457493837871696
ov2 = 0.016629818141012956
# semi_major_axis(or2,ov2)

tr2 = 1.054479678368431
tv2 = 0.016493452274361352
# semi_major_axis(tr2,tv2)

# semi_major_axis([or2,tr2],[ov2,tv2])


def eccentricity(v2,angular_mome,r2,mag_r2):
    if type(mag_r2) == list or type(mag_r2) == np.ndarray and len(mag_r2) == len(angular_mome):
        num = len(mag_r2)
    else:
        num = 1
    v2 = np.array(v2)
    angular_mome = np.array(angular_mome)
    r2 = np.array(r2)
    mag_r2 = np.array(mag_r2)
    mag_r2 = np.reshape(mag_r2,(num,1))
    e_vector = (np.cross(v2,angular_mome) / mue) - (r2/mag_r2)
    e = np.linalg.norm(e_vector,axis=1)  # unitless
    e = e.squeeze()
    e_vector = e_vector.squeeze()
    return e,e_vector

onv2 = [ 0.0032252372826048928 ,-0.01616118987443606 ,  -0.0006691444088957358]
oangm = [-0.0007609452114726807 ,-0.0008697921613171279  ,0.017339526789782626 ]
onr2 = [-1.0239826649099535 , -0.24517529866924212, -0.05723606361680865]
onmar2 = 1.054479678368431
# eccentricity(onv2,oangm,onr2,onmar2)

twv2 =[-0.0025254932789341334, -0.01629342444990578  , -0.002167268971870977 ]
tangm = [-0.003238414485647831 , -0.0021635186838188437,  0.020038916628412812 ]
twr2 = [-1.2298020637591551e+00 , 4.8688081059777844e-04, -1.9869115261357742e-01]
twmar2 = 1.2457493837871696
# eccentricity(twv2,tangm,twr2,twmar2)
# eccentricity([onv2,twv2],[oangm,tangm],[onr2,twr2],[onmar2,twmar2])

def inclination(angular_mo):
    angular_mo = np.array(angular_mo)
    if type(angular_mo[0]) == list or type(angular_mo[0]) == np.ndarray:
        angular_mo2 = angular_mo[:,2]
        num=1
    else:
        angular_mo2 = angular_mo[2]
        num=0
    mag_am = np.linalg.norm(angular_mo,axis=num)
    i = np.arccos(angular_mo2 / mag_am)
    i = np.degrees(i)
    i = i.squeeze()
    return i
an1 = [-0.0032384142967973456, -0.0021635186492247067,  0.020038916325615004 ]
# inclination(an1)
an2 = [-0.0007609454230707344 ,-0.0008697923613946281 , 0.017339526887897406 ]
# inclination(an2)
# inclination([an1,an2])

# functions assume an ecliptic frame of reference, so please if neccesary convert vectors before hand to get accurate results 
def longitude_ascending_node(angular_mo):
    angular_mo = np.array(angular_mo)
    z_hat = np.array([0,0,1])
    n = np.cross(z_hat,angular_mo)
    if type(angular_mo[0]) == list or type(angular_mo[0]) == np.ndarray and len(angular_mo[0]) == len(angular_mo[-1]):
        n0 = n[:,0]
        num = 1
        n1 = n[:,1]
    else:
        n0 = n[0]
        num = 0
        n1 = np.array(n[1])
    long_node = np.arccos((n0/np.linalg.norm(n,axis=num))) # can only return value from 0,pi
    long_node = np.array(long_node)
    long_node[n1 < 0] = 2 * np.pi - long_node[n1 < 0]
    long_node = np.degrees(long_node)
    return long_node

am1 = [-0.0032384142967973456, -0.0021635186492247067,  0.020038916325615004 ]
# longitude_ascending_node(am1)
am2 = [-0.0007609454230707344 ,-0.0008697923613946281 , 0.017339526887897406 ]
# longitude_ascending_node(am2)
# longitude_ascending_node([am1,am2])

def argument_periapsis(angular_mo,ecc):
    angular_mo = np.array(angular_mo)
    ecc = np.array(ecc)
    if type(angular_mo[0]) == list or type(angular_mo[0]) == np.ndarray and len(angular_mo[0]) == len(angular_mo[-1]):
        num = 1
        ecc2 = ecc[:,2]
    else:
        num = 0 
        ecc2 = ecc[2]
    z_hat = np.array([0,0,1])
    n = np.cross(z_hat,angular_mo)
    mag_ecc = np.linalg.norm(ecc,axis=num)
    mag_n = np.linalg.norm(n,axis = num)
    if type(angular_mo[0]) == list or type(angular_mo[0]) == np.ndarray and len(angular_mo[0]) == len(angular_mo[-1]):
        dotp = np.einsum('ij,ij->i', n, ecc)
    else:
        dotp = n@ecc
    
    arg_per = np.arccos(dotp/(mag_n*mag_ecc))
    arg_per = np.array(arg_per)
    arg_per[ecc2 < 0] = 2 * np.pi - arg_per[ecc2 < 0]
    arg_per = np.degrees(arg_per)
    return arg_per

mm1 = [-0.0007609454230707344, -0.0008697923613946281 , 0.017339526887897406 ]
ecc1 = [0.022116748634347094  ,0.045240265642180005  ,0.0032399543824206747]
# argument_periapsis(mm1,ecc1)
mm2 = [-0.0032384142967973456 ,-0.0021635186492247067 , 0.020038916325615004 ]
ecc2 = [-0.1320235265068338    , 0.1943515873911655    ,-0.0003524937083096147]
# argument_periapsis(mm2,ecc2)
# argument_periapsis([mm1,mm2],[ecc1,ecc2])

def mean_anomaly(mag_r,major_axis,ecc):
    mag_r = np.array(mag_r)
    major_axis = np.array(major_axis)
    ecc = np.array(ecc)
    ecc_anomaly = np.arccos((1-mag_r/major_axis)/ecc) 
    mean_anomaly =ecc_anomaly - (ecc * np.sin(ecc_anomaly))
    mean_anomaly = np.degrees(mean_anomaly)
    return mean_anomaly

mar1 = 1.245749369311783 
maj1 = 1.4905600329315598 
ec1 = 0.2349529215068008
# mean_anomaly(mar1,maj1,ec1)
mar2 = 1.0544796939365275 
maj2 = 1.0231605934888273
ec2 =0.050461168336950676
# mean_anomaly(mar2,maj2,ec2)

# mean_anomaly([mar1,mar2],[maj1,maj2],[ec1,ec2])

def orbital_elements(v2,r2):
    if type(v2[0]) == list or type(v2[0]) == np.ndarray and len(v2) == len(r2):
        v2 = np.array(v2,dtype=float)
        r2 = np.array(r2,dtype=float)
        num = 1
    else:
        v2 = np.array(v2)
        r2 = np.array(r2)
        num = 0
    r2 = rotate_from_equeritorial_to_ecliptic(r2)
    v2 = rotate_from_equeritorial_to_ecliptic(v2)
    mag_r2 = np.linalg.norm(r2,axis=num)
    mag_v2 = np.linalg.norm(v2,axis=num)
    angular_momentum = np.cross(r2,v2)
    a = semi_major_axis(mag_r2,mag_v2)# semi major axis
    e,e_vector = eccentricity(v2,angular_momentum,r2,mag_r2)
    i = inclination(angular_momentum)
    long = longitude_ascending_node(angular_momentum)
    peric = argument_periapsis(angular_momentum,e_vector)
    mean = mean_anomaly(mag_r2,a,e)
    return a,e,i,long,peric,mean

onev2 = [-0.002525492930422116, -0.014086836065775727, -0.008469579869150676]
oner2 = [-1.2298020505563299   ,0.07948145566433085 ,-0.18210191046123142]
# x = orbital_elements(onev2,oner2)
# print(x)
twov2 = [ 0.0032252366713229237 ,-0.014561432614395197 , -0.0070424777171607614]
twor2 = [-1.0239826861472883  ,-0.20217673754169574, -0.15003816245226842]
# x = orbital_elements(twov2,twor2)
# print(x)

# orbital_elements([onev2,twov2],[oner2,twor2])

def gauss_method(observation_times,p,R):
    """calculates the rhos and r values using the gauss method provided observatio_times, position and unit vector values 
    uses module math_functions,  and functions calculate_delta_t, calculate_a_and_b, calculate_cs, calculate_rhos_and_r

    Args:
        observation_times (list): list of modified julian date as [date1, date2, fate3]
        p (np.matrix): the unit matrix
        R (_type_): The position matrix

    Returns:
        list: a list with sublist that contains the caluclated roots, a list on weather the root is accepted or not, 
        all the r values for each of the roots, all the v2 values per each of the roots, and lastly all 
        the orbital elments where each column is a root and its rows the orbital elment. 
    """
    if type(observation_times[0]) == list or type(observation_times) == np.ndarray and len(p) == len(observation_times):
        size = len(observation_times)
    else:
        size = 1
    observation_times = np.array(observation_times)
    p = np.array(p)
    R = np.array(R)
    inverse_p = mf.take_inverse_matrix(p)
    B_matrix = inverse_p@R
    delta_t1,delta_t3,delta_t= calculate_delta_t(observation_times)
    a1,b1,a3,b3 = calculate_a_and_b(delta_t1,delta_t3,delta_t)
    d_1,d_2,magnitude_R_2,pos_2_dot_p2= calculate_values_for_cs(a1,b1,a3,b3,R,p,B_matrix)
    c,c3,c6,c8= calculate_cs(d_1,d_2,magnitude_R_2,pos_2_dot_p2)
    c8 = np.full(size, -1)
    c8 = c8.squeeze()
    roots = mf.eight_equation_four_coeff(c,c3,c6,c8) # need to fix what roots are removed...
    if size == 1:
        B_matrix = [B_matrix]
        a1 = [a1]
        b1 = [b1]
        a3 = [a3]
        b3 = [b3]
        R = [R]
        p = [p]
    new_B_matrix,new_a1,new_b1,new_a3,new_b3,new_R,new_p = [],[],[],[],[],[],[]
    for i in range(size):
        new_B_matrix.append([B_matrix[i].copy(),B_matrix[i].copy(),B_matrix[i].copy()])
        new_a1.append([a1[i].copy(),a1[i].copy(),a1[i].copy()])
        new_b1.append([b1[i].copy(),b1[i].copy(),b1[i].copy()])
        new_a3.append([a3[i].copy(),a3[i].copy(),a3[i].copy()])
        new_b3.append([b3[i].copy(),b3[i].copy(),b3[i].copy()])
        new_R.append([R[i].copy(),R[i].copy(),R[i].copy()])
        new_p.append([p[i].copy(),p[i].copy(),p[i].copy()])
    all_roots = []
    all_magnitudes = [] 
    all_r_values = []
    root_info = []
    for i in range(size):
        current_roots,magnitudes,r_values = calculate_rhos_and_r(roots[i],new_B_matrix[i],new_a1[i],new_b1[i],new_a3[i],new_b3[i],new_R[i],new_p[i])
        all_roots.append([current_roots])
        all_magnitudes.append(magnitudes)
        all_r_values.append(r_values)
        root_info.append([current_roots,magnitudes,r_values])
    all_roots = np.array(all_roots).squeeze()
    all_magnitudes = np.array(all_magnitudes).squeeze()
    all_r_values = np.array(all_r_values, dtype=np.matrix).squeeze()
    test_root = Check_root(root_info)
    test_root = np.array(test_root,dtype=bool).squeeze()
    v2 = []
    r2 = []
    all_orbit_elements = []
    all_r_values = all_r_values.squeeze()
    if size > 1:
        for j in range(size):
            v2_delta_t1 = np.repeat(delta_t1[j].copy(), 3)
            v2_delta_t3 = np.repeat(delta_t3[j].copy(), 3)
            current_v2 = calculated_v2(all_r_values[j],v2_delta_t1,v2_delta_t3)
            current_r2 = all_r_values[j,:,:,1]
            current_orbit_elements = (orbital_elements(current_v2,current_r2))
            r2.append(current_r2)
            all_orbit_elements.append(current_orbit_elements)
            v2.append(current_v2)
    else:
        v2_delta_t1 = np.repeat(delta_t1.copy(), 3)
        v2_delta_t3 = np.repeat(delta_t3.copy(), 3)
        v2 = calculated_v2(all_r_values,v2_delta_t1,v2_delta_t3)
        r2 = all_r_values[:,:,1]
        all_orbit_elements = orbital_elements(v2,r2)
        all_orbit_elements = np.array(all_orbit_elements)
    v2 = np.array(v2)
    r2 = np.array(r2)
    all_orbit_elements = np.array(all_orbit_elements)
    print(all_orbit_elements)
    print(all_roots)
    return all_roots,test_root,all_r_values,v2,all_orbit_elements
   
obs1 = [58577.48997074074, 58583.54555074074, 58590.54555074074]
p1 = [[-0.5258317338986549,  -0.5875255158277471  ,-0.6540863217081979 ],
      [ 0.8475382670837678 ,  0.8040126628785391  , 0.7481189653618031 ],
      [-0.07197133772397166 ,-0.09152817152276324, -0.11175463041961636]]
rr1 = [[-0.9696986007809081  ,-0.9406747492828259  ,-0.8945114818302049 ],
       [-0.22449591050121329, -0.316181379642978   ,-0.4177988254224926 ],
       [-0.09731285487753796 ,-0.13705996332878803, -0.18111530838131065]]
# gauss_method(obs1,p1,rr1)

obs2 = [58577.48997074074, 58601.54555074074, 58626.54555074074]
p2 = [[0.646743628403743,   0.5248759921239102 , 0.3627509124021436 ],
      [0.6000699731114835,  0.7060661103062899 , 0.8015212662825536 ],
      [0.47078520207111946, 0.4753691626187874 , 0.4753687360862347 ]]
rr2 = [[-0.9696986008561639,  -0.796540944160929   ,-0.4773994619640845 ],
       [-0.2244959102068138 , -0.564914570091051   ,-0.8191737144498188 ],
       [-0.09731285474991118 ,-0.24488569907572702, -0.3551090002088778 ]]
# gauss_method(obs2,p2,rr2)
# gauss_method([obs1,obs2],[p1,p2],[rr1,rr2])


def run_gauss_method_hand_values(date_time,obscode):
    """Uses function provided the values of the position and unit vecotr time and observetory code,
    it calculates the gauss_method and orbital elements uses function provide_unit_matrix and 
    provide_position_matrix (they are expected to alreydy have the corect matrices)
    
    Args:
        date_time (list): a list of the times
        obscode (string): observatory code

    Returns:
        list: a list of maximum of lenght three (the lenght depends on the number of real positive roots generated), 
        where is made up of sublist with the following order
        [[root, [[rho1],[rho2],[rho3]], [[1x3 matrix r1],[1x3 matrix r2],[1x3 matrix r3]]], [root2,...]...]
    """
    postion_matrix = provide_position_matrix() 
    unit_matrix = provide_unit_matrix() 
    result = gauss_method(date_time,unit_matrix,postion_matrix)
    return result

# run_gauss_method_hand_values([58577.489970740738,58583.545550740739,58590.545550740739],'500')

def run_gauss_method_wis(date_time,ra_dec_values,obscode):
    """Uses function calculate_position_obs, get_unit_matrix, gauss_method
    
    Args:
        date_time (list): a list of the times
        ra_dec_values (list): list where each elements is a sublist with ra and dec values, in the follwoing format:
        [[Ra1,Dec1],[Ra2,Dec2]...] 
        obscode (string): observatory code

    Returns:
        list: a list of maximum of lenght three (the lenght depends on the number of real positive roots generated), 
        where is made up of sublist with the following order
    """
    postion_matrix = calculate_position_obs(date_time.copy(),obscode) # seperated from gauss, to make calculation faster
    postion_matrix = postion_matrix.T
    unit_matrix = get_unit_matrix(ra_dec_values) # [[alpha,delta][]...]
    result = gauss_method(date_time,unit_matrix,postion_matrix)
    return result

# run_gauss_method_wis([58577.489970740738,58583.545550740739,58590.545550740739],[[2.1260971174823160,-7.2033616739254860E-002],[2.2018577151544165,-9.1656450482163324E-002],[2.2892342062190609,-0.11198856664053504]],'500')

# asteroid 17188   
ra_dec31 = [0.74798128583882861,0.49018056855941849]
ra_dec32 = [0.93154498092039317,0.49538358898508605]
ra_dec33 = [1.1458005817864867,0.49538310417140496]
all_three_val = [ra_dec31,ra_dec32,ra_dec33]
time2 = [58577.489970740738,58601.545550740739,58626.545550740739]
# run_gauss_method_wis(time2,all_three_val,'500')

test21_ra_dec = [1.1536801615828238,0.20914571314856717]
test22_ra_dec = [1.2043302709741688,0.21933455747076516]
test23_ra_dec = [1.2623223260301843,0.22963781782170506]
test2_times = [58577.489970740738,58583.545550740739,58590.545550740739]
# run_gauss_method_wis(test2_times,[test21_ra_dec,test22_ra_dec,test23_ra_dec],'500')
