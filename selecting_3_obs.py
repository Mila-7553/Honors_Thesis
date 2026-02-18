# print("hello world")
import numpy as np
import random
import gauss.gauss_method as gm

# Gathering dates
test_dates = [1.8492,2.9382,3.8402,4.3828,5.8482,6.8482,6.8382,7.032,7.329,8.5654,8.53,9.1,9.542,10.23,282,28371,1737,27291,2.32,1.423,3443,12]
print(len(test_dates))
exit()
# Sorting dates
# doing something with a lot of dates given
def create_A_3obs_sets(dates: np.array, total_arc_lenght=15,max_mid_last = 5,max_mid_first=5,min_mid=2):
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
    date_obs.sort() 
    sets_3_obs = []
    for i in range(len(date_obs)):    
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
                        sets_3_obs.append(current_selection)
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
            print("this is subspace")
            print(subspaces[-1])
            start = end
        # select a random values from each of the subspaces
        selected_values = np.array([random.choice(subspace) for subspace in subspaces])
    else:
        # When indexes are directly provided by the user, those indesex are choosen to retun the sets of dates. 
        try:
            selected_values = np.array(obs_sets[indexes])
        except:
            print("Please provide valid indexes for choosing sets of dates")
            print("selecting_3_obs.py, chose_obs_sets")
            raise ValueError
    return selected_values

x = create_A_3obs_sets(test_dates)

sets = chose_obs_sets(x,5)
print(len(sets))
print(sets)

# next step, get real dates and testing, + making sure WIS integration 
# can work with sets of 3 observations that belong to different observatories