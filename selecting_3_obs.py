# print("hello world")
import numpy as np
import random
import gauss.gauss_method as gm
from astropy.time import Time

test_dates = [1.8492,2.9382,3.8402,4.3828,5.8482,6.8482,6.8382,7.032,7.329,8.5654,8.53,9.1,9.542,10.23,282,28371,1737,27291,2.32,1.423,3443,12]

def create_A_3obs_sets(dates: np.array, total_arc_lenght=15,max_mid_last = 5,max_mid_first=5,min_mid=2,max_comb=None):
    """Makes a list of sets of 3 observations, with all possible combinations that are 
    accepted based on the prvided constraints passed as parameters

    Args:
        dates (np.array): a list of dates to be sorted, and generated valid sets of 3 dates
        total_arc_lenght (int, optional): The maximum arc lenght that is accepted. Defaults to 15.
        max_mid_last (int, optional): maximum (time) lenght between the second observation and 
        the third observation. Defaults to 5.
        max_mid_first (int, optional): maximum (time) lenght between the second observation
        and the first observation. Defaults to 5.
        min_mid (int, optional): minimal distance that the first and last observation 
        must be from the second observation. Defaults to 2.
        max_comb (None/int): provides the maximum number of combinations that can be return/ Defaults to None

    Raises:
        ValueError: if the provided dates, are not valid floats 

    Returns:
        np.array: array containning all sets of 3 of all valid combinations of the provided dates
    """
    try:
        date_obs = np.array(dates, dtype=np.float64)
    except:
        print("Please provide valid dates for making sets of 3 observations")
        print("Selecting_3_obs.py, create_A_3obs_sets")
        raise ValueError
    # chekning that there are less then 1000 obs provided
    if len(date_obs) >= 1000:
        print(f"Too many observations provided {len(date_obs)}")
        raise ValueError
    date_obs.sort() 
    sets_3_obs = []
    count = 0 # total number of sets of three found
    for i in range(len(date_obs)):    
        if float(date_obs[i]) <= 2451545.5000000:
            continue
        for j in range(i+1, len(date_obs)):
            diff_obs2_obs1 = date_obs[j] - date_obs[i]
            
            if (diff_obs2_obs1 > max_mid_first):
                break
            elif (diff_obs2_obs1 >= min_mid):
                
                for k in range(j+1, len(date_obs)):
                    diff_obs3_obs1 = date_obs[k] - date_obs[i]
                    diff_obs3_obs2 = date_obs[k] - date_obs[j]
                    
                    if (diff_obs3_obs1 > total_arc_lenght) or (diff_obs3_obs2 > max_mid_last):
                        break
                    elif (diff_obs3_obs2 >= min_mid):
                        current_selection = [date_obs[i], date_obs[j], date_obs[k]]
                        count += 1
                        # using reservoir sampling
                        # putting a limit in the man number of combination stores based on passed argument
                        if max_comb is None: # when no limit on number of possible combination
                            sets_3_obs.append(current_selection)

                        else:
                            if len(sets_3_obs) < max_comb:
                                sets_3_obs.append(current_selection)
                            else:
                                r = np.random.randint(0, count)
                                if r < max_comb:
                                    sets_3_obs[r] = current_selection
                        # gives all combiantion equal probablilty to be part of the selected sets of three
    # worst case O(n^3), wll be slow with very loose constraints 
    sets_3_obs = np.array(sets_3_obs, dtype=np.float64)
    return sets_3_obs

def chose_obs_sets(obs_sets: np.array,num_sets=10,indexes=[]):
    """provide with a list of multiple sets of 3 observations randomly and 
    symetrically chooses the specifide number of sets. 

    Args:
        obs_sets (np.array): contains all sets of 3 observations to choose from
        num_sets (int, optional): How many set you want it to return. Defaults to 10.
        indexes (list, optional): if you know the indexes of the sets that you want to
        choose you can manually provide them. Defaults to [].

    Raises:
        ValueError: when provided indexes is not a valid number (index out of range or not a number)

    Returns:
        np.array: radomly selcted sets from the provided list of sets. 
    """
    if not indexes:
        # select at random, from similarly spaced spaces/sets, (subspace)
        total_sets = len(obs_sets)
        subspace_size = total_sets // num_sets # decide normal size of a subspace
        remainder = total_sets % num_sets # number of additions that will be done to first couple of subspaces
        # 12 element per subspace
        # 4 reminder
        # Making subspace with the same as the number of sets
        subspaces = []
        start = 0 # initiation first subspace on the first index
        for i in range(num_sets):
            # Adjust subspace size for remainder
            j = 0 # for keeping track of if one (index) should be added to the
            # current subspace (if there is a reminder)
            if remainder != 0:
                j = 1
                remainder -=1
            end = start + subspace_size + j
            subspaces.append(obs_sets[start:end])
            start = end
        # select a random values from each of the subspaces
        selected_values = []
        for quad in subspaces:
            value = random.choice(quad)
            # inn = np.where(obs_sets == value)[0]
            selected_values.append(value)
        selected_values = np.array(selected_values)
    else:
        # When indexes are directly provided by the user, those indesex are choosen to retun the sets of dates. 
        try:
            selected_values = np.array(obs_sets[indexes])
        except:
            print("Please provide valid indexes for choosing sets of dates")
            print("selecting_3_obs.py, chose_obs_sets")
            raise ValueError
    return selected_values

# print(create_A_3obs_sets(test_dates))
# x = create_A_3obs_sets(test_dates,max_comb=10)

# print(x)
# exit()

# sets = chose_obs_sets(x,5)
# print(len(sets))
# print(sets)

# next step, get real dates and testing, + making sure WIS integration 
# can work with sets of 3 observations that belong to different observatories

def position_obs_wis(mjd,obscode): # uses Wis.py
    """Uses Wis library to calculate the position of the observatory in a Heliocentric equatorial frame. 
    Provided the mjd and observatory code. the mjd can be provided as sets of dates, however the 
    observatory code can only be a single string

    Args:
        mjd (float or list): list of dates of the observation 
        obscode (string): Observatory code of the observation

    Returns:
        np.matrix: a 3x3 matrix per each of the given times 
        (such as provided a list of 4 times, it returns a list with 4, 3x3 matrices corresponding to each time)
    """
    # flattening the parameters 
    jd_times = gm.flat_variable(mjd)
    obscode = gm.flat_variable(obscode)
    jd_times = gm.mjd_to_jd(jd_times)
    
    if type(obscode) == str or type(obscode) == int or len(obscode) == 1:
        obscode = np.repeat(obscode,len(mjd))
    
    if len(obscode) != len(jd_times):
        raise ValueError("Please provide the same number of observatory codes and times")
    # grouping observatory codes and corresponding times in a dictonary
    dict_code_times = {}
    for i in range(len(obscode)):
        if obscode[i] not in dict_code_times:
            dict_code_times[obscode[i]] = []
        dict_code_times[obscode[i]].append([jd_times[i],i])  # putting index on the dictonary for reodering 
   
    positions = np.zeros((len(mjd),3))
    from wis import wis  as wi
    for code, values in dict_code_times.items():
        values = np.array(values)
        jd_times = values[:,0]
        indexes = np.array(values[:,1],dtype=int)
        times = Time(jd_times, format='jd', scale='tdb')
         # imported withing the function, because  importing takes some runtime
        W = wi(str(code), times) # as w and pass the w
        observer_posns = W.hXYZ
        observer_posns = np.array(observer_posns)
        positions[indexes] = observer_posns
    return positions
    # return np.matrix(observer_posns)
obscodes = ["A24","703","505","500","500", "505","500"]
times = [58577.489970740738,58583.545550740739,58590.545550740739]
times = [58577.489970740738,58577.489970740738,58577.489970740738,58577.489970740738,58577.489970740738,58577.489970740738,58577.489970740738]
# print(len(times))
# print()

# observation has, time, asteroid, ra, dec, obscode 

# Class Observation

# What it intially returns [3x3 matrix]
# [[-0.9696986008561639  -0.2244959102068138  -0.09731285474991118] onw row corresponds to one time and one obs code
#  [-0.9406747493915482  -0.31618137935660046 -0.13705996320464234]
#  [-0.8945114819710074  -0.4177988251583791  -0.18111530826682148]]
# x= position_obs_wis(times,obscodes)
# print(x)

# list1 = np.zeros((10,3))
# indexes = np.array([1,5,2])
# values = np.matrix([[43,543,21],[324,439,439],[321,232,565]])
# list1[indexes]=values
# print(list1)
# print()
# print(x)
# 17188        KC2012 12 11.18706 05 23 22.28 +33 56 15.7          18.1 Rc~0nALA24
# 17188         C2012 12 20.22671 05 00 48.16 +30 43 54.6          17.2 Vr~0mv3703
