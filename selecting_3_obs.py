print("hello world")
import numpy as np
# Gathering dates
dates = [1.8492,2.9382,3.8402,4.3828,5.8482,6.8482,6.8382,7.032,7.329,8.5654,8.53,9.1,9.542,10.23,282,28371,1737,27291,2.32,1.423,3443,12]

# Sorting dates
date_obs = np.array(dates, dtype=np.float64)
date_obs.sort()

max_constraint_arc_length = 15 # five days constrain on total arc length
max_constranit_mid_last = 3
max_constraint_mid_first = 10 
min_constraint_mid = 1 # The smallets accept range from the the mid obs
# slow 3 loops no cutting aproach 
sets_3_obs = []
for i in range(len(date_obs)):
    
    for j in range(i+1, len(date_obs)):
        diff_obs2_obs1 = date_obs[j] - date_obs[i]
        
        if (diff_obs2_obs1 > max_constraint_mid_first):
            break
        elif (diff_obs2_obs1 >= min_constraint_mid):
            
            for k in range(j+1, len(date_obs)):
                diff_obs3_obs1 = date_obs[k] - date_obs[i]
                diff_obs3_obs2 = date_obs[k] - date_obs[j]
                
                if (diff_obs3_obs1 > max_constraint_arc_length) or (diff_obs3_obs2 > max_constranit_mid_last):
                    break
                elif (diff_obs3_obs2 >= min_constraint_mid):
                    current_selection = [date_obs[i], date_obs[j], date_obs[k]]
                    sets_3_obs.append(current_selection)
# worst case O(n^3), wll be slow with very loose constraints 
sets_3_obs = np.array(sets_3_obs, dtype=np.float64)