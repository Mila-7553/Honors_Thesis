import numpy as np
from gauss import gauss_method as gm
import pyorb
from astropy.time import Time, TimeDelta
import astropy.units as u
import astropy.constants as const
from scipy.optimize import least_squares

# calculation the orbital elments. 
eclip = np.radians(23.439281)# used for coordinate rotation
def rotate_equirotarial_to_ecliptic(vect):
    """This functions converts the given vector to ecliptic coordindates. 
    the function can perform the calculations for multiple object if parameters are given as list 

    Args:
        vect (list): a list of at least lenght three made up of vectors of lenght
        three that are going to be rotates to ecliptic coordinate

    Returns:
        list: The vector that was given in the parameter but in ecliptic coordinates
    """
    vect = np.array(vect)
    vect_total = vect.flatten()
    vect_total = int(len(vect_total) / 3) # determining the amount of sets
    
    # making sure the dimension are okay
    if vect.ndim > 2:
        vect = np.reshape(vect,(vect_total,3))
    
    if vect.ndim == 1:
        vect = np.reshape(vect, (1,3))
    
    # math    
    rot_matrix = np.array([[1,0,0],[0,np.cos(-eclip),-1*np.sin(-eclip)],[0,np.sin(-eclip),np.cos(-eclip)]])
    rotation = np.array([np.copy(rot_matrix) for _ in range(vect_total)])
    vect = np.reshape(vect, (vect_total, 3,1))
    res = rotation @ vect
    
    # cleaning up solution
    res = np.reshape(res, (vect_total, 3))
    res = res.squeeze()
    return res

def semi_major_axis(r2,v2):
    """Function used to calculate the semi major axis an object, this function does 
    not correctly calculate the semi major asixis for object with hyperbolic orbits.
    it can perform the calculations for multiple object if parameters are given as list

    Args:
        v2 (_type_): v2 vector of the object, as an array of lengh 3
        angular_mome (_type_): angular momentum of the object
        r2 (list): r2 vector of the object, as an arry of lengh 3
        mag_r2 (float): magnitude of the r2 vector of the object
    Returns:
        list / float: returns the semi major axis for the object(s) provided
    """
    # flattening the provided parameter, and cheking their sizes conscide
    r2 = gm.flat_variable(r2)
    v2 = gm.flat_variable(v2)
    if len(r2) != len(v2):
        print("Error, semi_major_ais function on the file gauss_method.py")
        print("please, check the dimension of the values that you provided")
        print(len(r2),len(v2))
        raise ValueError
    
    # math
    a = ((2/r2)-(v2**2/gm.mue))**(-1) # units AU
    a = a.squeeze()
    return a

def eccentricity(v2,angular_mome,r2,mag_r2):
    """Function used to calculate the eccentricity and the eccentricity vector of an object,
    it can perform the calculations for multiple object if parameters are given as list
    
    Args:
        v2 (_type_): v2 vector of the object, as an array of lengh 3
        angular_mome (_type_): angular momentum of the object
        r2 (list): r2 vector of the object, as an arry of lengh 3
        mag_r2 (float): magnitude of the r2 vector of the object

    Returns:
        tuple: Returns the eccetricity and the eccentricity vector for the object(s) provided
    """
    # Flattens and convers the variables into numpy arrays
    mag_r2 = gm.flat_variable(mag_r2)
    size = len(mag_r2)
    v2 = np.array([v2])
    angular_mome = np.array([angular_mome])
    r2 = np.array([r2])
    
    # changes the dimension of the parameters v2 and r2 if neccesary
    if v2.ndim > 2:
        v2 = np.reshape(v2,(size,3))
        angular_mome = np.reshape(angular_mome,(size,3))
        r2 = np.reshape(r2,(size,3))
        
    # math
    mag_r2 = np.array(mag_r2)
    mag_r2 = np.reshape(mag_r2,(size,1))
    e_vector = (np.cross(v2,angular_mome) / gm.mue) - (r2/mag_r2)
    e = np.linalg.norm(e_vector,axis=1)  # unitless
    e = e.squeeze()
    e_vector = e_vector.squeeze()
    return e,e_vector

def inclination(angular_mo):
    """Function used to calculate the inclination an object, it can perform
    the calculations for multiple object if parameters are given as list
    
    Args:
        angular_mo (list): The angular mommentum of the object as a vector of lenght three

    Returns:
        float / list: Returns the argument periapsis for the object(s) provided
    """
    # Flattens the provided parameter and determines the amount of sets.
    angular_mo = np.array(angular_mo)
    size = angular_mo.flatten()
    size = int(len(size) / 3)
    
    # changes the dimensions of the angular momnetum if neccesary
    if angular_mo.ndim > 2 or angular_mo.ndim < 2:
        angular_mo = np.reshape(angular_mo,(size,3))
    
    #math
    angular_mo2 = angular_mo[:,2]
    mag_am = np.linalg.norm(angular_mo,axis=1)
    i = np.arccos(angular_mo2 / mag_am)
    i = np.degrees(i)
    i = i.squeeze()
    return i

def longitude_ascending_node(angular_mo):
    """Function used to calculate the longitude ascending node of an object, it can perform
    the calculations for multiple object if parameters are given as list
    
    Args:
        angular_mo (list): The angular mommentum of the object as a vector of lenght three

    Returns:
        float / list: Returns the longitude ascending node for the object(s) provided
    """
    # flattens the given parameter, and detemines the number of sets present as size
    angular_mo = np.array(angular_mo)
    size = gm.flat_variable(angular_mo)
    size = int(len(size) / 3)
    
    # changes the dimension of the angular momentum parameter if necessary
    if angular_mo.ndim > 2 or angular_mo.ndim < 2:
        angular_mo = np.reshape(angular_mo,(size,3))
        
    # math
    z_hat = np.array([0,0,1])
    n = np.cross(z_hat,angular_mo)
    n0 = n[:,0]
    n1 = n[:,1]
    long_node = np.arccos((n0/np.linalg.norm(n,axis=1))) # can only return value from 0,pi
    long_node = np.array(long_node)
    long_node[n1 < 0] = 2 * np.pi - long_node[n1 < 0] 
    long_node = np.degrees(long_node)
    long_node = long_node.squeeze()
    return long_node

def argument_periapsis(angular_mo,ecc_vector):
    """Function used to calculate the argument periapsis of an object, it can perform
    the calculations for multiple object if parameters are given as list

    Args:
        angular_mo (list): The angular mommentum of the object as a vector of lenght three
        ecc_vector (list): The eccentricity vector of lenght 3 

    Returns:
        float / list : Returns the argument periapsis for the object(s) provided
    """
    
    # makes the provided parameters numpy arrays and flattens them to check their 
    # sizes concsides 
    ecc = np.array(ecc_vector)
    angular_mo = np.array(angular_mo)
    size_ecc = gm.flat_variable(ecc)
    size_angm = gm.flat_variable(angular_mo)
    size_ecc = int(len(size_ecc) / 3)
    size_angm = int(len(size_angm) / 3)
    if size_angm != size_ecc:
        print("Error, argument_periapsis function on the file gauss_method.py")
        print("please, check the dimension of the values that you provided")
        print(size_ecc,size_angm)
        print(ecc.shape,angular_mo.shape)
        raise ValueError
    
    # Cheking the dimension of the angular momnetum and changing the dimension if neccesary.
    if angular_mo.ndim > 2 or angular_mo.ndim < 2:
        angular_mo = np.reshape(angular_mo,(size_angm,3))
        ecc = np.reshape(ecc,(size_angm,3))
    ecc2 = ecc[:,2]
    
    # math
    z_hat = np.array([0,0,1])
    n = np.cross(z_hat,angular_mo)
    mag_ecc = np.linalg.norm(ecc,axis=1)
    mag_n = np.linalg.norm(n,axis = 1)
    dotp = np.einsum('ij,ij->i', n, ecc)
    arg_per = np.arccos(dotp/(mag_n*mag_ecc))
    arg_per = np.array(arg_per)
    arg_per[ecc2 < 0] = 2 * np.pi - arg_per[ecc2 < 0]
    arg_per = np.degrees(arg_per)
    arg_per = arg_per.squeeze()
    return arg_per

def mean_anomaly(mag_r,major_axis,ecc):
    """Function used to calculate the mean anamoly of an object, it can perform
    the calculations for multiple object if parameters are given as list

    Args:
        mag_r (list / float): The magnitude of the r vector provied.
        major_axis (list / float): the major axis of the object
        ecc (list / float): the eccentricity of the object

    Returns:
        np.ndarray: the mean anomoly of the provided objects
    """
    # flattening parameters and cheking their size consciede
    mag_r = gm.flat_variable(mag_r)
    major_axis = gm.flat_variable(major_axis)
    ecc = gm.flat_variable(ecc)
    if len(mag_r) != len(major_axis):
        print("Error, mean_anomaly function on the file gauss_method.py")
        print("please, check the dimension of the values that you provided")
        print(len(mag_r), len(major_axis), len(ecc))
        raise ValueError
    
    # math
    ecc_anomaly = np.arccos((1-mag_r/major_axis)/ecc) 
    mean_anomaly =ecc_anomaly - (ecc * np.sin(ecc_anomaly))
    mean_anomaly = np.degrees(mean_anomaly)
    mean_anomaly = mean_anomaly.squeeze() # removing any additionall brackets
    return mean_anomaly, ecc_anomaly

def orbital_elements(v2,r2,eclip=0):
    """Given vectors of length 3, `v2` and `r2`, this function calculates the orbital elements  
    of an asteroid: semi-major axis (a), eccentricity (e), inclination (i), longitude of the  
    ascending node (long), argument of periapsis (peric), and mean anomaly (mean).  

    It can calculate the orbital elements for either a single set of vectors or multiple sets.  

    Args:  
        v2 (list): A list (length divisible by 3) representing one or more v2 vectors.  
        r2 (list): A list (length divisible by 3) representing one or more r2 vectors.
        eclip: if it is set to zero rotates the provided vrctors from equirotarial 
        to ecliptic coordinates to pereform calculations.
        
    Returns:
        list: Orbital elements for each set of provided v2 and r2 values, structured 
            as separate arrays for each element type in the order [a, e, i, long, peric, mean].  
            Each array contains the values for all sets in order.  
            Example: [[a1, a2, a3], [e1, e2, e3], [i1, i2, i3], ...].
    """
    # making the parameters flatten numpy arrays, and cheking their sizes conscide 
    v2 = np.array(v2)
    r2 = np.array(r2)
    size_v2 = gm.flat_variable(v2)
    size_r2 = gm.flat_variable(r2)
    if len(size_r2) != len(size_v2):
        print("Error, orbital_elements function on the file gauss_method.py")
        print("please, check the dimension of the values that you provided")
        print(v2.shape, r2.shape)
        print(len(size_v2),len(size_r2))
        raise ValueError
        
    size = int(len(size_r2) /3) # number of sets
    
    # changing the coordinates to perform math efficiently
    if eclip == 0:
        r2 = rotate_equirotarial_to_ecliptic(r2)
        v2 = rotate_equirotarial_to_ecliptic(v2)
    
    # cheking and changing the dimension of the parameters if neccesary
    if v2.ndim > 2 or v2.ndim < 2:
        v2 = np.reshape(v2,(size,3))
        r2 = np.reshape(r2,(size,3))
        
    # math and function calls to calculate orbital elements
    mag_r2 = np.linalg.norm(r2,axis=1)
    mag_v2 = np.linalg.norm(v2,axis=1)
    angular_momentum = np.cross(r2,v2)
    a = semi_major_axis(mag_r2,mag_v2)
    e,e_vector = eccentricity(v2,angular_momentum,r2,mag_r2)
    i = inclination(angular_momentum)
    long = longitude_ascending_node(angular_momentum)
    angular_momentum = angular_momentum.squeeze()
    peric = argument_periapsis(angular_momentum,e_vector)
    mean,ecc_mean = mean_anomaly(mag_r2,a,e)
    
    # extra feature for when eclip != used for testing of functionalitty 
    # and scientific accuracy.
    if eclip != 0:
        mean = ecc_mean
        peric = np.radians(peric)
        long = np.radians(long)
        i = np.radians(i)    
    return a,e,i,long,peric,mean


def unit_calc_ra_dec(unit_vector):
    # needs input as [x,y,z],[x,y,z],...
    unit_vector = np.array(unit_vector)
    # print(unit_vector)
    try:
        mag_unit = np.linalg.norm(unit_vector,axis=1)
    except np.exceptions.AxisError:
        mag_unit = np.linalg.norm(unit_vector)
        unit_vector = np.array([unit_vector])

    # print("is this unit vector: ", unit_vector)
    ra = np.arctan2(unit_vector[:,1],unit_vector[:,0]) * 180 / np.pi
    ra = gm.flat_variable(ra)
    # atan2(vEQ[1],vEQ[0])*180./PI 
    dec = (np.pi / 2 - np.acos(unit_vector[:,2]/mag_unit))*180 / np.pi
    dec = gm.flat_variable(dec)
    for i in range(len(ra)):
        if ra[i] < 0:
            ra[i] += 360
        if ra[i] > 360:
            ra[i] -= 360
        if dec[i] < -90:
            ra[i] = dec[i] + 180
        if dec[i] > 90:
            ra[i] = dec[i] - 180
    ra = np.radians(ra)
    dec = np.radians(dec)
    return np.array([ra,dec])

# differential correction Still on progress
def make_pyorb_orbit(r_vec, v_vec, epoch, G):
    return pyorb.Orbit(
        M0=pyorb.M_sol,
        epoch=epoch,
        G=G,
        x=r_vec[0],
        y=r_vec[1],
        z=r_vec[2],
        vx=v_vec[0],
        vy=v_vec[1],
        vz=v_vec[2]
    )

def my_propagate(dt,orbit):
        if not isinstance(dt,TimeDelta):
            raise TypeError(f"{type(dt)=}: But it is supposed to be a TimeDelta object")

        # Update epoch as an astropy Time
        # orbit.epoch = orbit.epoch + dt

        # Call parent propagate using dt in its expected units (here: days)
        return orbit.propagate(dt.to_value(u.day))

def get_cartesian_at_times(orbit, og_epoch, times):
    N = len(times)
    XYZ = np.zeros((N,3)) # initializing array that will contain xyz components

    # making a copy of the orbit
    orb = orbit.copy()
    time = np.squeeze(times)
    for i in range(N):
        # exit()
        dt = TimeDelta(time[i] - orb.epoch)  # TimeDelta
        # neworb = my_propagate(dt,orb)
        orbit.epoch = orbit.epoch + dt
        my_propagate(dt,orb)
        XYZ[i,:] = orb.r[:,0]
        # print(neworb.r)
        # exit()
    return XYZ

    
def light_time_correction(orbit, epoch, obs_times, obs_positions):
    N = len(obs_times)
    obs_positions = np.asarray(obs_positions)

    # starting with no light delay
    delay = np.zeros(N)
    obs_times = gm.flat_variable(obs_times)
    i = 0 # counting the # of iterations
    while True:
        i += 1
        # time of emission = time of observation - delay
        
        emit_times = (obs_times - delay)
        # position of asteroid at emission time
        obj_pos = get_cartesian_at_times(orbit, epoch, emit_times)

        # observer → object vectors
        sepn = obj_pos - obs_positions
        dist = np.linalg.norm(sepn, axis=1)

        # new delay (AU divided by AU/day = days)
        new_delay = dist / (const.c.to('AU/day').value)
        
        # check convergence
        if np.allclose(new_delay, delay) or i > 200:
            break

        delay = new_delay
     # unit vectors
    return sepn / dist[:,None]
    
def ecliptic_to_equatorial(vecs):
    """
    Rotate ecliptic vectors to equatorial J2000.
    vecs: Nx3
    """
    eps = np.deg2rad(23.439291111)   # obliquity of ecliptic J2000
    ce = np.cos(eps)
    se = np.sin(eps)

    R = np.array([
        [1, 0, 0],
        [0, ce, -se],
        [0, se,  ce]
    ])
    return vecs @ R.T
    
def get_ra_dec_lightTime_correction(r2, v2, obs_times, obs_positions, G=gm.mue):

    # making the epoch the second time of observation
    epoch = obs_times[1]

    orbit = make_pyorb_orbit(r2, v2, epoch, G)

    ecl_vectors = light_time_correction(
        orbit, epoch,
        obs_times,
        obs_positions
    )
    e = unit_calc_ra_dec(ecl_vectors)
    print('Ra and Dec value with no conversion: ',e)
    eq_vectors = ecliptic_to_equatorial(ecl_vectors)
    # print(eq_vectors)
    x = unit_calc_ra_dec(eq_vectors)
    print('Ra and Dec value after conversion from ecliptic to equatorial: ',x)

    # return ra_dec


def get_ra_dec_no_lightTime_correction(r_vectors, obs_positions):
    
    ra_dec = []
    unit_vectors = []
    for r, obs in zip(r_vectors, obs_positions):
        rho = np.array(r) - np.array(obs) # [r_ix - obs_ix, r_iy - obs_iy, r_iz - obs_ix] 
        # lenght of unit vectro from center of the sun - 
        # lenght of position of the observaroty = leght rho 
        unit_vect_rho = rho / np.linalg.norm(rho) # direction of rho as a unit vector
        # conversion
        unit_vectors.append(unit_vect_rho)
        ra = np.arctan2(unit_vect_rho[1], unit_vect_rho[0]) % (2*np.pi)
        dec = np.arcsin(unit_vect_rho[2])
        ra_dec.append([ra, dec])
    return np.array(ra_dec), unit_vectors


def residual(predicted_unit_vector,starting_unit_vector,observer_position):
    observer_position = observer_position.T
    predicted_unit_vector = np.array(predicted_unit_vector)
    if predicted_unit_vector.ndim == 1:
        predicted_unit_vector = np.reshape(predicted_unit_vector,(3,3))
    
    resid = np.degrees(np.arccos(np.clip(np.sum(starting_unit_vector * predicted_unit_vector, axis=-1), -1.0, 1.0)))
    print("this is residual: ", resid,flush=True)
    return resid

result = least_squares(
        residual,  # Function to compute residuals
        np.array(unit).T.flatten(),  # Initial guess, flattened
        args=(knonw_unit_vector,position),  # Arguments passed to the residual function
        method='trf',  # Optimization method
        verbose=2  # Verbosity level
    )
