# import pyorb
# from astropy.time import Time, TimeDelta
# import astropy.units as unit
from gauss import gauss_method as gm
import numpy as np
import pyorb
from astropy.time import Time, TimeDelta
import astropy.units as u
import astropy.constants as const
# import healpy

# ============================================================
# (0) BUILD ORBIT FROM CARTESIAN + EPOCH
# ============================================================

def make_orbit_cartesian(r_vec, v_vec, epoch, G):
    """
    Create a pyorb orbit from Cartesian state at a known epoch.
    r_vec, v_vec are 3-component numpy arrays (AU, AU/day).
    epoch is an astropy Time object.
    """
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
    

# ============================================================
# (1) PROPAGATE TO A GROUP OF TIMES
# ============================================================

def get_cartesian_at_times(orbit, epoch, times):
    """
    Compute Cartesian r(t) for multiple times.
    orbit  – pyorb.Orbit (mutable!)
    epoch  – astropy.Time, the epoch of r2,v2 (middle observation)
    times  – astropy.Time array (times of observation)
    """
    N = len(times)
    XYZ = np.zeros((N,3))

    # We must not mutate the original orbit → make a copy
    orb = orbit.copy()
    time = np.squeeze(times)
    for i in range(N):
        dt = TimeDelta(time[i] - epoch)  # TimeDelta
        orb.epoch = epoch            # reset epoch each time
        orb.propagate(dt.to_value(u.day))
        XYZ[i,:] = orb.r[:,0]

    return XYZ


# ============================================================
# (2) LIGHT-TIME ITERATION
# ============================================================

def light_time_corrected_vectors(orbit, epoch, obs_times, obs_positions):
    """
    Computes the apparent direction from observer to asteroid,
    including light-time correction.

    obs_times: length-N astropy Time array
    obs_positions: Nx3 array of observatory Cartesian position (AU)
    """
    N = len(obs_times)
    obs_positions = np.asarray(obs_positions)

    # initial guess: no light delay
    delay = np.zeros(N)
    obs_times = gm.flat_variable(obs_times)
    i = 0
    while True:
        i += 1
        # time of emission = time of observation - delay
        
        emit_times = (obs_times - delay) * u.day
        # position of asteroid at emission time
        obj_pos = get_cartesian_at_times(orbit, epoch*u.day, emit_times)

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


# ============================================================
# (3) CONVERT ECLIPTIC UNIT VECTORS → EQUATORIAL → RA, DEC
# ============================================================

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


# def vec_to_radec(equatorial_vectors):
#     """
#     Convert Nx3 equatorial vectors to RA,Dec (radians).
#     Uses healpy.vec2ang.
#     healpy returns angles as (theta,phi) or (lon,lat).
#     We request lonlat=True → (RA_deg, Dec_deg)
#     """
#     radec_deg = np.array(healpy.vec2ang(equatorial_vectors, lonlat=True)).T
#     return np.deg2rad(radec_deg)     # convert to radians


# ============================================================
# (4) COMPLETE PIPELINE
# ============================================================

def cartesian_to_radec(vecs):
    """
    Convert equatorial Cartesian vectors to RA, Dec.
    vecs: shape (N,3)
    Returns RA,Dec in radians.
    """
    vecs = np.array(vecs)

    if vecs.ndim == 1:
        vecs = vecs.reshape(1,3)

    x = vecs[:,0]
    y = vecs[:,1]
    z = vecs[:,2]

    r = np.sqrt(x*x + y*y + z*z)

    # Declination
    dec = np.arcsin(z / r)

    # Right ascension
    ra = np.arctan2(y, x)
    ra = np.mod(ra, 2*np.pi)

    return np.vstack((ra, dec)).T

def compute_RaDec_from_roots(r_vectors, obs_positions):
    """
    Compute RA/Dec directly from FORTRAN-style Gauss roots.
    r_vectors: 3x3 array, object vectors (r1, r2, r3)
    obs_positions: 3x3 array, observatory positions
    """
    ra_dec = []
    for r, obs in zip(r_vectors, obs_positions):
        vec = np.array(r) - np.array(obs)
        unit_vec = vec / np.linalg.norm(vec)
        ra = np.arctan2(unit_vec[1], unit_vec[0]) % (2*np.pi)
        dec = np.arcsin(unit_vec[2])
        ra_dec.append([ra, dec])
    return np.array(ra_dec)

import numpy as np

import numpy as np

def compute_RaDec_from_roots_lighttime_iterative(r_vectors, obs_positions, times, r_dot2, t2, c=299792.458, tol=1e-8, max_iter=10):
    """
    Compute RA/Dec including iterative light-time correction.
    
    r_vectors: 3x3 array, object vectors at observation times (initial guess)
    obs_positions: 3x3 array, observatory positions
    times: observation times corresponding to r_vectors
    r_dot2: velocity at middle time t2
    t2: middle time
    c: speed of light (km/s)
    tol: convergence tolerance in distance
    max_iter: maximum iterations for light-time convergence
    """
    ra_dec = []
    r2 = r_vectors[1]  # middle position
    
    for r_obs, obs, t_obs in zip(r_vectors, obs_positions, times):
        # initial guess: observed position
        r_corr = np.array(r_obs)
        
        for _ in range(max_iter):
            # vector from observer to object
            vec = r_corr - obs
            delta_t = np.linalg.norm(vec) / c
            # propagate from middle position back to emission time
            r_new = r2 + r_dot2 * (t_obs - t2 - delta_t)
            # check convergence
            if np.linalg.norm(r_new - r_corr) < tol:
                break
            r_corr = r_new
        
        unit_vec = (r_corr - obs) / np.linalg.norm(r_corr - obs)
        ra = np.arctan2(unit_vec[1], unit_vec[0]) % (2*np.pi)
        dec = np.arcsin(unit_vec[2])
        ra_dec.append([ra, dec])
    
    return np.array(ra_dec)



# def compute_RaDec_from_gauss(r2, v2, obs_times, obs_positions, G=gm.mue):
#     """
#     Compute RA/Dec for three observations using Gauss orbit.
    
#     r2, v2         – Cartesian state at middle observation
#     obs_times      – 3-element astropy.Time array (t1, t2, t3)
#     obs_positions  – 3x3 array of observer positions (AU)
#     G              – gravitational parameter
#     """
#     epoch = obs_times[1]  # middle observation

#     # Build orbit at middle observation
#     orbit = make_orbit_cartesian(r2, v2, epoch, G)

#     # Propagate to all observation times
#     r_at_times = get_cartesian_at_times(orbit, epoch, obs_times)

#     # Compute unit vectors and RA/Dec
#     ra_dec = []
#     for r_obj, obs in zip(r_at_times, obs_positions):
#         vec = r_obj - np.array(obs)          # vector from observer to object
#         unit_vec = vec / np.linalg.norm(vec)
#         ra = np.arctan2(unit_vec[1], unit_vec[0]) % (2*np.pi)
#         dec = np.arcsin(unit_vec[2])
#         ra_dec.append([ra, dec])

#     ra_dec = np.array(ra_dec)
#     print("RA/Dec in radians:\n", ra_dec)
#     print("RA/Dec in degrees:\n", np.rad2deg(ra_dec))

#     return ra_dec

# # ============================================================
# # (4) EXAMPLE INPUT
# # ============================================================

position =[ 
    [-0.96969860078090808, -0.22449591050121329, -0.09731285487753796],
    [-0.84476964478826122, -0.50018798005564036, -0.21682977412483659],
    [-0.61876282670557059, -0.73286918395739420, -0.31769850609934919]
]

r2 = [-0.83606191165353250, -0.48940308600177118, -0.20935702854324587]
r_vector =[ [-0.96025912405227454    ,  -0.21573765150910806 ,      -9.0441567436148257E-002],
           [-0.83606191165353250  ,    -0.48940308600177118  ,    -0.20935702854324587],
           [-0.61221659642978876 ,     -0.72126273475966607   ,   -0.31047690349056062]]
v2 = [8.9775543621979960E-003, -1.3255504452997285E-002, -5.7691400089880072E-003]

times = [
    gm.mjd_to_jd(58577.489970740738),
    gm.mjd_to_jd(58596.545550740739),
    gm.mjd_to_jd(58616.545550740739)
]

ra = compute_RaDec_from_roots(r_vector, position)
print(ra)
print("here")
# exit()
def compute_RaDec_from_gauss(r2, v2, obs_times, obs_positions, G=gm.mue):
    """
    r2, v2         – Cartesian state found by Gauss’s method (middle observation)
    obs_times      – length-3 astropy.Time array (t1,t2,t3)
    obs_positions  – 3x3 array of observatory position vectors (AU)
    G              – gravitational parameter used by pyorb

    RETURNS:
        unit vectors (ecliptic, light-time corrected)
        RA, Dec in radians
    """

    # --------------------------------------------------------
    # Epoch = time of middle observation  ← THIS IS CRITICAL
    # --------------------------------------------------------
    epoch = obs_times[1]

    # Build the orbit at epoch
    orbit = make_orbit_cartesian(r2, v2, epoch, G)

    # Light-time corrected vectors (in ecliptic frame)
    ecl_vectors = light_time_corrected_vectors(
        orbit, epoch,
        obs_times,
        obs_positions
    )

    # Convert to equatorial frame
    eq_vectors = ecliptic_to_equatorial(ecl_vectors)
    ra_dec = cartesian_to_radec(eq_vectors)
    print("this is ra dec resulting with propagation:", ra_dec)

    return ecl_vectors, eq_vectors

position =[ [-0.96969860078090808  ,    -0.22449591050121329     ,  -9.7312854877537963E-002],
           [-0.84476964478826122  ,    -0.50018798005564036     , -0.21682977412483659],
           [-0.61876282670557059  ,    -0.73286918395739420,      -0.31769850609934919 ]]
r2 = [-0.83606191165353250     , -0.48940308600177118  ,    -0.20935702854324587]
v2 = [8.9775543621979960E-003 , -1.3255504452997285E-002,  -5.7691400089880072E-003]
times = [gm.mjd_to_jd(58577.489970740738),gm.mjd_to_jd(58596.545550740739),gm.mjd_to_jd(58616.545550740739)]
print()
print('here')
x = compute_RaDec_from_gauss(r2,v2,times,position)
print(x)
print()
exit()
print("_"*10)
import numpy as np

def compute_RaDec_iterative_lighttime(r2, r_dot2, t2, obs_positions, obs_times, c=299792.458, tol=1e-12, max_iter=50):
    """
    Compute RA/Dec including iterative light-time correction.

    Parameters
    ----------
    r2 : np.ndarray
        3-element array of position at middle time (reference) [same units as obs_positions, e.g. km]
    r_dot2 : np.ndarray
        3-element array of velocity at middle time [same units per unit time]
    t2 : float
        Middle observation time (same units as obs_times)
    obs_positions : np.ndarray
        Array of shape (N,3) with observer positions at observation times
    obs_times : np.ndarray
        Array of N observation times corresponding to obs_positions
    c : float
        Speed of light in same units as r2 per unit time (default km/s)
    tol : float
        Convergence tolerance for light-time (distance units)
    max_iter : int
        Maximum number of iterations for light-time convergence
    """
    ra_dec_list = []

    for obs_pos, t_obs in zip(obs_positions, obs_times):
        # initial guess: emission position = middle position
        r_emission = r2.copy()

        for _ in range(max_iter):
            vec = np.array(r_emission) - np.array(obs_pos)
            dist = np.linalg.norm(vec)
            delta_t = dist / c
            r_new = r2 + r_dot2 * (t_obs - t2 - delta_t)
            if np.linalg.norm(r_new - r_emission) < tol:
                break
            r_emission = r_new

        # Compute unit vector from observer to object
        vec_unit = (r_emission - obs_pos) / np.linalg.norm(r_emission - obs_pos)

        # RA / Dec
        ra = np.arctan2(vec_unit[1], vec_unit[0]) % (2 * np.pi)
        dec = np.arcsin(vec_unit[2])

        ra_dec_list.append([ra, dec])

    return np.array(ra_dec_list)

ee = compute_RaDec_iterative_lighttime(r2,v2,times[1],position,times)
print(ee)
e = compute_RaDec_from_roots_lighttime_iterative(r_vector,position,times,v2,times[1])
print(e)
# exit()
vector =[ 0.23031989828824631    ,   0.87871604357837652      , 0.38094652534937923]
# Edit r2_1:   0.20320988178720450     
# Edit r2_2:    1.5777860402453072     
# Edit r2_3:   0.80651462837170340
vector_r2 = [0.20320988178720450 ,1.5777860402453072,0.80651462837170340]
vector_v2 = [-1.0598295382117201E-002 , -4.1509093492139246E-003 , -8.1529719286995497E-003]
# r = gm.unit_calc_ra_dec(vector)
vector_times = [56269.352817592589,56281.227487592594,56284.192173592593]
# print(r)
print("_"*10)
# compute_RaDec_from_gauss(vector_v2,vector_v2,vector_times,vector)
exit()

import numpy as np

def ra_deg_format(ra_rad):
    """
    Convert RA (radians) into degrees and arcminutes:
    e.g. 53.28367° → '53° 17.02′'
    """
    # convert to degrees
    ra_deg = np.rad2deg(ra_rad)
    
    # normalize 0–360
    ra_deg = ra_deg % 360.0

    d = int(ra_deg)
    arcmin = (ra_deg - d) * 60.0

    return f"{d:02d}° {arcmin:05.2f}′"


def dec_deg_format(dec_rad):
    """
    Convert Dec (radians) into degrees, arcminutes, arcseconds:
    e.g. 29.4753° → '+29° 28′ 31.1″'
    """
    # convert to degrees
    dec_deg = np.rad2deg(dec_rad)

    sign = "+" if dec_deg >= 0 else "-"
    dec_abs = abs(dec_deg)

    d = int(dec_abs)
    m_float = (dec_abs - d) * 60.0
    m = int(m_float)
    s = (m_float - m) * 60.0

    return f"{sign}{d:02d}° {m:02d}′ {s:04.1f}″"


def radec_format(ra_rad, dec_rad):
    """Return formatted RA and Dec strings."""
    return ra_deg_format(ra_rad), dec_deg_format(dec_rad)

ra_dec = [1.314453277730429, 0.39708404892679416]
# 1.4397776446635013     
# Edit alpha 2:   1.3124994704239334     
# Edit alpha 3:   1.2794529728748587     
# Edit delta 1:   0.60626716568484362     
# Edit delta 2:   0.53637163239656471
ra_dec = [1.3124994704239334,0.53637163239656471]
x = radec_format(ra_dec[0],ra_dec[1])
print(x)
#  jd_time = gm.mjd_to_jd(58583.545550740739)
# obj_orbit = make_orbit_obj(-1.2298020505563314,0.07948145566433279,-0.18210191046123164,-0.0025254929304212755,-0.014086836065777309,-0.008469579869150574,Time(jd_time,format='jd'))

# def make_orbit_obj(x, y, z, vx, vy, vz,epoch, G=gm.mue):
#     # epoch is the time of the observation
#     orbit = pyorb.Orbit(M0=pyorb.M_sol, epoch=epoch, G=G, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
#     return orbit

# def propogate(orbit, epoch, dt):
#     """_summary_

#     Args:
#         orbit (_type_): _description_
#         epoch (_type_): _description_
#         dt (_type_): _description_

#     Returns:
#         _type_: _description_
#     """
#     if not type(dt) == TimeDelta:
#         raise TypeError(f"{type(dt)=}: But it is supposed to be a TimeDelta object")
#     # orb = pyorb.Orbit()
#     # orb.from_elements(orbit[0], orbit[1], orbit[2], orbit[3], orbit[4], orbit[5], epoch)
#     # new_orb = orb.propagate(dt)
#     # epoch= epoch + dt
#     orbit.epoch = epoch + dt
#     orbit.propagate(dt.to_value(unit.day))

# jd_time = gm.mjd_to_jd(58583.545550740739)
# obj_orbit = make_orbit_obj(-1.2298020505563314,0.07948145566433279,-0.18210191046123164,-0.0025254929304212755,-0.014086836065777309,-0.008469579869150574,Time(jd_time,format='jd'))

# propogate(obj_orbit, Time(jd_time,format='jd'), TimeDelta(7.0, format='jd'))


# def light_po

# def create_orbit_state(M0, epoch: Time, G, x, y, z, vx, vy, vz, **kwargs):
#     """
#     Creates a new OrbitState-like structure without using a class.
#     M0, G, x, y, z, vx, vy, vz are parameters for the orbital state.
#     """
#     if not isinstance(epoch, Time):
#         raise TypeError(f"'epoch' must be an astropy.time.Time instance, got {type(epoch)}")

#     if not epoch.isscalar:
#         raise ValueError(f"'epoch' must be a scalar Time (single value), got length {len(epoch)}")

#     orbit = Orbit(M0=M0, epoch=epoch, G=G, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, **kwargs)
#     return orbit, epoch

# def propagate(orbit, epoch: Time, dt: TimeDelta):
#     """
#     Propagate orbit over a time-span dt (in days or as a Quantity convertible to days),
#     and advance epoch consistently.
#     """
#     if not isinstance(dt, TimeDelta):
#         raise TypeError(f"{type(dt)=}: But it is supposed to be a TimeDelta object")

#     # Update epoch as an astropy Time
#     epoch = epoch + dt

#     # Call pyorb.Orbit's propagate method using dt in its expected units (here: days)
#     orbit.propagate(dt.to_value(u.day))
#     return orbit, epoch

# def get_XYZs(orbit, epoch: Time, times: TimeDelta):
#     """
#     Propagate the orbit to the supplied times, and extract the Cartesian positions.
#     """
#     if times.isscalar:
#         times = Time([times])  # Make iterable if necessary

#     XYZs = np.empty((len(times), 3), dtype=float)
#     for i, time in enumerate(times):
#         orbit, epoch = propagate(orbit, epoch, time - epoch)
#         XYZs[i, :] = orbit.r.T[0]
#     return XYZs

# def generate_unit_vector(orbit, epoch: Time, times: Time, observatoryXYZ: np.ndarray, approx=False):
#     """
#     Calculate apparent UnitVector from specified observatory-posn(s) to orbit-object at given time(s).
#     """
#     if not isinstance(times, Time):
#         raise TypeError("'times' must be an astropy.time.Time instance")

#     if times.isscalar:
#         times = Time([times])  # Make iterable if necessary

#     # Ensure that the shapes of the input data are consistent
#     if observatoryXYZ.shape != (len(times), 3):
#         raise ValueError(f"{observatoryXYZ.shape=} != {(len(times), 3)=}")

#     # Initialize light-delay to 0
#     lightDelay = np.zeros(len(times))
#     while True:
#         # Calculate delayed time (of emission)
#         objectXYZ = get_XYZs(orbit, epoch, times - lightDelay * u.day)

#         # Calculate relative separation-vector from observatory to object
#         sepn_vectors = objectXYZ - observatoryXYZ

#         # Calculate distance to object at each time
#         d = np.linalg.norm(sepn_vectors, axis=1)

#         # Calculate light-travel-time
#         thisLightDelay = d / (astropy.constants.c * 86400 / astropy.constants.au).value

#         # Break out of loop if converged
#         if np.allclose(thisLightDelay, lightDelay) or approx:
#             break
#         lightDelay = thisLightDelay

#     # Return unit-vectors: shape = (N_times, 3)
#     return sepn_vectors / d[:, None]

# def ecliptic_to_equatorial(ecliptic_vectors):
#     """
#     Rotate ecliptic vectors to equatorial vectors
#     """
#     if ecliptic_vectors.ndim != 2 or ecliptic_vectors.shape[1] != 3:
#         raise ValueError("Input ecliptic_vectors must have shape (N, 3)")

#     ecl = (84381.448 * (1. / 3600) * np.pi / 180.)  # Obliquity of ecliptic at J2000
#     ce = np.cos(ecl)
#     se = np.sin(ecl)
#     rotmat = np.array([[1.0, 0.0, 0.0],
#                        [0.0, ce, -se],
#                        [0.0, se, ce]])
#     return rotmat.dot(ecliptic_vectors.T).T

# def ecliptic_UV_to_RaDec(ecliptic_vectors):
#     """
#     Convert ecliptic vectors to (RA, Dec) in an equatorial coordinate frame
#     """
#     equatorial_vectors = ecliptic_to_equatorial(ecliptic_vectors)
#     return np.radians(np.array(healpy.vec2ang(equatorial_vectors, lonlat=True)).T)

# def new_from_cartesian(cartesian, M0, epoch, G):
#     """
#     Create new orbit state from cartesian coordinates.
#     """
#     if cartesian.ndim != 1 or cartesian.shape[0] != 6:
#         raise ValueError(f'new_from_cartesian: Expected shape (6,), got {cartesian.shape=}')
#     orbit, epoch = create_orbit_state(M0, epoch, G, *cartesian)
#     return orbit, epoch
