# Honors_Thesis
Gauss method for orbital determination code.
## download the code:
To download the code clone the repository in your local machine, within the directory where you cloned your repository, go to the terminal to the directory of your repository ...\Honors_Thesis and run the command  *python3 -m pip install .*
additionally make sure the following python libraries are installed: 
- astropy 
- bs4             
- numpy            
- pandas           
- pip               
- pytest             
- scipy             
## install and use the code:
Once the code has download and all the requirements are set, go to the folder gauss open the file fortran_code_replica and run it, this will run a demonstration of the code, and If it does not crash the code was successfully install in your machine, if it does crash please make sure the code was properly downloaded. 
## run the tests:
For each of the files within the gauss directory there exit a corresponding test file inside the test directory, to access and use this file, in the terminal go to the directory containing the github repository test folder ...\Honors_Thesis\test, once in this directory run the command *python3 -m pytest*  This runs all the test on all the files of the directory at once. To run the test file in a single file use *python3 -m pytest "file_name.py"*. 
Additionally Github automatically runs all the test and on each pull and push, and you can also manually trigger it to run the test, the test results can be seen under actions. 