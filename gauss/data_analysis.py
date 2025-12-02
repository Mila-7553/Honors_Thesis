import sqlite3
from pathlib import Path
from gauss import gauss_method as gm
from gauss import fortran_code_replica as fr
from gauss import additional_features as adf
import numpy as np

def get_cart_kep(lines):
    """provided a set of lines searches for the numerical values from lines containing
    'CAR'(cartesian coordinates) and 'KEP'(keplerian elements) and returns the results as lists

    Args:
        lines (list):a list of string to check for lines containg 'CAR' and 'KEP'

    Returns:
        tuple: containes two lists: of the keplerian elements and cartesian values.
    """
    cart_values = []
    kep_values = []
    for line in lines:
        if 'CAR' in line:
            my_lines= line.strip().split(" ")
            for i in my_lines:
                try:
                    num = float(i)
                    cart_values.append(num)
                except:
                    pass 
            if len(cart_values) == 12: 
                cart_values = cart_values[:6]
                
        if line.strip().startswith('KEP'):
            my_lines = line.strip().split(' ')
            for i in my_lines:
                try:
                    num = float(i)
                    kep_values.append(num)
                except:
                    pass 
            if len(kep_values) == 12: 
                kep_values = kep_values[:6]
    return kep_values, cart_values

def read_fel_ke(name):
    """This function reads a file .ke from orbfit resulting output files, and returns 
    the keplerian and cartesian data of the file, uses function get_cart_kep.

    Args:
        name (str): The name of the .ke file (shuuld include the .ke ending)

    Returns:
        tuple: a tuple where the frst element is the keplerian data found on the file 
        and the second element in the cartesian coordinates foind on the file
    """ 
    script_dir = Path(__file__).resolve().parent # getting parent directory
    file_path = script_dir / "fortran_data" / name 
    with open(file_path, "r") as f:
        my_lines = f.readlines()
    kep, cart = get_cart_kep(my_lines)
    return kep, cart

def read_calc(name):
    """This function reads a file containing data and extracts the roots from the line 
    containign accpeterd root values from each line to return a list of accpeted roots.

    Args:
        name (string): name of the file

    Returns:
        list: list of accepted roots on the file
    """
    script_dir = Path(__file__).resolve().parent
    file_path = script_dir / "fortran_data" / name
    with open(file_path, "r") as f:
        roots = []
        for line in f:
            if 'accepted root' in line:
                my_lines= line.strip().split(" ")
                current_root = []
                for i in my_lines:
                    try:
                        num = float(i)
                        current_root.append(num)
                    except:
                        pass 
                roots.append(current_root[:2])
    return roots

# functions create_db, edit_db, make_table and insert_info, 
# came from previusly written code by me in class COMP-390 software engineering
def create_db(name):
    """
    Creates a new file for the database with the name provided, it erases all data previusly store in that file
    if it alredy existed, and set ups a connection with the database
    :param name: a string, that is the name of the file storing the database
    :return: a tuple first element is the connection to the database and the second element is the cursor
    """
    # making sure database file is empty so data is not place twice if the program is run more then once
    with open(name, "w") as db:
        db.close()
    # establishing the connection with the database
    conn = sqlite3.connect(name)
    cursor = conn.cursor()
    return conn, cursor

def edit_db(name):
    """
    Creates a conection to an existing database.
    
    :param name: a string, the name of the file storing the database
    :return: a tuple where the first element is the connection to the database 
    and the second element is the cursor.
    """
    try:
        conn = sqlite3.connect(name)
        cursor = conn.cursor()
    except sqlite3.Error as e:
        print(f"Error accesing database: {e}")
        return None

    return conn, cursor

def make_table(table_name, colums_specification, cursor, conn):
    """
    If a table does not exist in the database, it will create it, with columns provided by the user in the lis
    colums_specification. This function also checks if we were able to make the table successfully if this is
    not the case it prints an error messages AND continues

     :param table_name: The name of the table being created
     :param colums_specification: a list of the colums that you want to use with all specification included such as
     reclong TEXT (name of the colum and which data type shoul be used)
    """
    colums = ""
    for i in range(len(colums_specification)):
        if i != len(colums_specification) - 1:
            colums += colums_specification[i] + ", "
        else:
            colums += colums_specification[i]

    sql = f"""CREATE TABLE IF NOT EXISTS {table_name} ({table_name}_ID INTEGER PRIMARY KEY AUTOINCREMENT,{colums});"""
    try:
        # try and except stament for any time we attempt to connect to the database
        cursor.execute(sql)
        conn.commit()
    except:
        # error message that will exit the program if the function fails
        print(f"We were unable to create the table {table_name}")
        print(
            f"please make sure, that the table_name is in the correct format, "
            f"and that the database is connected correctly"
        )
        print(sql)
        exit()

def insert_info(colums_data, table_name, colums, cursor, conn):
    """
    provided a dictonary with its keys representing colum names with an EXISTING table, the values of the dictonary
     would be the correspoinding value that gets added to the databses such as {'colum1':'data1'} there is a table with
    a colum named colum1, for which we will star a new value ina row below colum1, we place data 1
    :param colums_data: a dictonary where the keys are colum names, and values of the dictonary, respective values
    pased to the database
    :param table_name: a string containg the name of the table that we want to insert the data to
    :return: None
    """
    # not possible to have or check for nulls (all provided list by the user should, have all nulls included
    names = colums
    # putting the names into string format that sql understands
    # for key in colums_data.keys():
    #     # Not all dictonary have the same keys so we do a check
    #     if key in colums:
    #         names.append(key)
    #     else:
    #         remove.append(key)
    # colums_data = remove_from_dict(
    #     remove, colums_data
    # )  # removing information from the dictonary
    
    colum_names = "','".join(names)  # making it a string
    colum_names = f"'{colum_names}'"  # adding ' at the start and end of the string
    values = ()
    for (val) in (colums_data.values()):  # a .join was not used to unite the data since if it was use the values
        if type(val) == str:
            # if the value is a string it checks that characters that could cause issue with getting the data to the
            # sql, and replaces them with a similar looking charcter
            val = val.replace("'", "`")
            val = val.replace('"', "`")
        if type(val) == dict or type(val) == list:
            # sql query gets confuse when provided a list or dictonary of data so it was converted to a string
            val = str(val)
        values += (val,)
    paramet = ""
    for i in range(len(colums_data.keys())):
        # since only one row is added at a time, I realize that each parameter needs a ?
        # this this for loop creates a ? for each value
        paramet += "?,"
    paramet = paramet[:-1]  # removes extra comma at the end of the list of ?,...,
    # values = (values[3])
    sql = (
        f"INSERT INTO {table_name} " f"({colum_names}) " f"VALUES ({paramet});"
    )  # making the complete insert statement
    # makes sure that the sql is succefully executed, if it is not, it prints a messages AND continues the program

    try:
        cursor.execute(sql, values)
        conn.commit()
    except:
        print(f"We were unable to insert the data into the table {table_name}")
        print(f"please make sure, that the data is in an appropriate format data: ")
        print(f"make sure that the table {table_name} exist in the database")
        print("and that the database is correctly connected")
        print(sql)
        print(values)
        exit()
        return
    return

def load_data_db():
    """This function uses functions create db, to make a connection with a new db,
    function make_table, to create on the connected db, functions read_calc, read_fel_ke, 
    to read file and extracr specific useful information from those files to 
    then use function insert_info topopulates a database table with the data extracted
    """
    conn, cursor =create_db("test_db.db",)
    colums_specification = [ 'Asteroid VARCHAR(200)','Asteroid_root INT', 'Correction INT', 'a INT', 'e INT', 'i INT', 'long_node INT', 'peric INT', 'mean_anomaly INT', 'x INT', 'y INT', 'z INT', 'vx INT', 'vy INT', 'vz INT', 'comments VARCHAR(200)']
    make_table("fortran_data", colums_specification, cursor, conn)
    files_names = ['s17188','s00433','153814_tracklet_35_40','17188_tracklet_35_40']
    for file_name in files_names:
        roots = read_calc(f'Calc_values_{file_name}.txt')
        for num, root in roots:
            ast_name = file_name.split('_')[0]
            if ast_name.startswith('s'):
                ast_name = file_name[1:]
            # o means no correction, 1 means with correction
            dict_asteroid = {'Asteroid':f'{ast_name}', 'Asteroid_root':root, 'Correction':'0','comments':f'file name: {file_name}.obs, root no. {num}, no correction'}
            kep_val, cart_val = read_fel_ke(f'r{int(num)}_{file_name}.ke')
            labe_cart = ['x','y','z','vx','vy','vz']
            label_kep = ['a','e','i','long_node','peric','mean_anomaly']
            ke_dict_kep = dict(zip(label_kep, kep_val))
            ke_dict_cart = dict(zip(labe_cart, cart_val))
            ke_merged_dict = dict_asteroid | ke_dict_kep | ke_dict_cart
            insert_info(ke_merged_dict, "fortran_data", ke_merged_dict.keys(), cursor, conn)
            kep_val, cart_val = read_fel_ke(f'r{int(num)}_{file_name}.fel')
            # exit()
            dict_asteroid['Correction'] = 1
            dict_asteroid['comments'] = f'file name: {file_name}.obs, root no. {num}, with correction'
            labe_cart = ['x','y','z','vx','vy','vz']
            label_kep = ['a','e','i','long_node','peric','mean_anomaly']
            fel_dict_kep = dict(zip(label_kep, kep_val))
            fel_dict_cart = dict(zip(labe_cart, cart_val))
            fel_merged_dict = dict_asteroid | fel_dict_kep | fel_dict_cart
            insert_info(fel_merged_dict, "fortran_data", fel_merged_dict.keys(), cursor, conn)
    conn.commit()
    conn.close()

def file_line_chek_basic(line):
    """This function parses different types of information from a given line and returns a numerical
    identifier and corresponding value, uses function.

    Args:
        line (String): line of text as input and extracts certain values based on 
        specific keywords present in the line.

    Returns:
        list: returns a tuple containing two values: `num` and `value`. 
        The `num` variable represents a numerical code indicating the type of information extracted
        from the input line, while the `value` variable holds the extracted value based on the content of
        the input line.
    """
    value = ''
    num = 0 
    if 'alpha' in line:
        num = 1
        value = np.float64(line.strip().split(" ")[-1])
        
    if 'delta' in line:
        num = 2
        this_delta = line.strip().split(" ")[-1]
        value = np.float64(this_delta)
        
    if 'time' in line:
        num = 3
        this_time = line.strip().split(" ")[-1]
        value = np.float64(this_time)
        
    if 'Position' in line:
        num = 4
        current_vect = np.array(line.strip().split(" "))
        value_indexes = current_vect[:] != ''
        current_vect = current_vect[value_indexes]
        value = np.array(current_vect[-3:],dtype=np.float64)
    if "Edit root" in line:
        num = 5
        value = np.float64(line.strip().split(" ")[-1])
        
    if "Edit r:" in line:
        num = 6
        r_values = np.array(line.strip().split(" "))
        value_indexes = r_values != ''
        r_values = r_values[value_indexes]
        value = np.array(r_values[2:], dtype=np.float64)
        
    if "Edit vp" in line:
        num = 7
        vp = np.array(line.strip().split(" "))
        value_indexes = vp != ''
        vp = vp[value_indexes]
        value = np.array(vp[2:], dtype=np.float64)
        
    if 'accepted' in line:
        num = 100
        
    if 'spurious root' in line:
        num = -1
        
    return num, value 

def organize_line_values(flag, vp, ra, r, roots, obsP, times,dec):
    """This function organizes line values based on certain conditions and returns the organized
    information, without including the flag parameter.

    Args:
        flag (lis): indicates true if the root [i] index is a spurious root or false otherwise
        vp (list): a list of velocity of the asteroid at r_2 when the roots is not 
        spurious containg [vp1, vp2,vp3] where the lenght of the list varies depending 
        on if the root is purous or no and where vpn is [vpn_x, vpn_y, vpn_z]
        ra (list): List of right ascencion 
        r (list): list of the vector positions at r1, r2 and r3, as a [r1,r2,r3] where 
        rn is represented by [rn_x,rn_y, rn_z]
        roots (list): a list of roots 
        obsP (list): list of positions as provided by fortran 
        [position1, postion2, position3] where each postion_n is = [xn,yn,zn]
        times (list): list of times
        dec (list): list of declenation values

    Returns:
        list: a list in the following order [ra,dec,time,observerposition,root1,root2,root3] 
        where each rootn contains [root_value,r_values of the root,vp if no spurious]
    """
    info = [ra,dec,times,obsP]
    for i in range(len(flag)):
        my_root_info = []
        my_root_info.append(roots[i])
        my_root_info.append(r[i])
        if flag[i] == False:
            current_vp = vp.pop(0)
            my_root_info.append(current_vp)
        info.append(my_root_info)
    return info

def check_line_by_line(lines):
    """this function uses functions file_line_chek_basic and organize_line_values,
    to take an input of lines, and get the root information and other information from this lines

    Args:
        lines (list): The lines that are going to be checked

    Returns:
        list: the values extracted from the lines, in the follwoing 
        format [ra,dec,time,observerposition,root1,root2,root3] 
        where each rootn contains [root_value,r_values of the root,vp if no spurious]
    """
    spurious_flag, vp_values, ra_values, r_values, roots, obs_positions, times,dec_values = ([] for _ in range(8))
    for line in lines:
        my_num,value = file_line_chek_basic(line)
        if my_num == -1:
            spurious_flag.append(True)
            # we have a spurius root therefore values for vp do not get score
        elif my_num == 100:
            spurious_flag.append(False)
        if my_num == 1:
            ra_values.append(value)
        if my_num == 2:
            dec_values.append(value)
        if my_num == 3:
            times.append(value)
        if my_num == 4:
            obs_positions.append(value)
        if my_num == 5:
            roots.append(value)
        if my_num == 6:
            r_values.append(value)
        if my_num == 7:
            vp_values.append(value)
    my_values = organize_line_values(spurious_flag, vp_values, ra_values, r_values, roots, obs_positions, times,dec_values)        
    return my_values

def get_initial_data(file_name):
    """This function reads a file located in a subfolder named "fortran_data" and then
    checks its contents line by line using function check_line_by_line.

    Args:
        file_name (str): the name of the file to be checked
    
    Returns:
        list: the result from using check_line_by_line function, the values extracted from the lines,
        in the follwoing format [ra,dec,time,observerposition,root1,root2,root3] 
        where each rootn contains [root_value,r_values of the root,vp if root is not spurious]
    """
    script_dir = Path(__file__).resolve().parent
    file_path = script_dir / "fortran_data" / file_name
    with open(file_path, "r") as f:
        file_lines = f.readlines()
    my_data = check_line_by_line(file_lines)
    return my_data
    

files = ['calc_values_17188_tracklet_1_12.txt','calc_values_17188_tracklet_35_40.txt','calc_values_153814_tracklet_35_40.txt','calc_values_s00433.txt','calc_values_s17188.txt']
for i in files:
    i = 'calc_values_s00433.txt'
    file_info = get_initial_data(i)
    current_ra = file_info[0]
    current_dec = file_info[1]
    current_obs_times = file_info[2]
    current_obs_pos = np.array(file_info[3]).T
    root1 = file_info[4]
    
    current_unit_vector = gm.unit_vector_from_ra_dec(current_ra,current_dec)
    x = gm.gauss_method(current_obs_times,current_unit_vector,current_obs_pos)
    orb = adf.orbital_elements(x[-2][0],x[-1][0])
    print(orb)
    exit()
    

