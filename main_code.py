#-------------------------------------------------------------------------------------------------------------------
# This file is the main code from where all the other files are accessed. To start the whole optimisation process and 
# the pipeline, only this file need to be run from Python, and the rest will be accessed as desired.
#
# The pipeline assumes a Mooney-Rivlin three parameter material model is used, with the C01, C10 and C20 as the
# parameters. Also it is assumed that 3D displacement field is used.

# The code does account for 2D and also the Mooney-Rivlin two parameter model with parameters C01 and C10. But might 
# contain some minor bugs and is not guaranteed to run smoothly.

# The whole pipeline can also be used to test two optimisation algorithms. During such a procedure only the type of
# algorithm need to be changed in line 273 . Therefore, this file will run for the 1st algorithm, there after the 
# a few changes are needed in the lines and files below and then this file need to be run again:
#       -   "main_code.py"
#               - lines 179 - 202 explain what needs to be done in this file.
#       - "graphs.py"
#       - "File_paths.py"
#
# Francè Bresler
# 12 February 2020
# -------------------------------------------------------------------------------------------------------------------

# All the packages and modules needed for the file
import math as m
import numpy as nm
import dot11 as dot
from py_mentat import *
from py_post import *
import os
import subprocess
import pandas as pd
import csv
import sys
import time
from functions import material2d, material3d
from File_paths import filepaths
from code_settings import pipeline_files
pf = pipeline_files()
# ------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------
# The "myEvaluate" function, is what DOT use to evaluate the objective function and contraints.
# The RMS error is calculated and used as the objective function
# A boolean contraint is used to tell DOT if inequality constraint was adhered to or not.
#   - if the MARC analysis converged, g[0] = -1
#   - if the MARC analysis didn't converge, g[0] = 1
def myEvaluate(x, obj, g, param):
    #Select which Marc Mentat version to use for the analysis
    mentat_version = filepaths(pf.mentat_version)
    # Since 10 deisgn points will be optimised, the optimser need to know which design point's 
    # optimisation is currently running, the "iterations" file store the current design point's 
    # number, zero indexed.
    # The current deisgn point's number only gets accessed here.
    file = open("iterations.txt","r")
    it = int(file.readline())
    file.close
    # //

    # //
    # A bias exits between the 3 design variables, therefore all of them gets scaled to a value
    # of 1 by dividing the starting design point by itself. This eliminates the bias within the 
    # optimisation calculations.
    # For each optimisation iteration the design variables need to be multiplied again by the 
    # starting point to ensure the procedure file written to MARC, runs the FE analysis with the 
    # correct variables. 
    # The next few lines obtaines the starting point to be multiplied with the current design variables.
    path1 = filepaths("path1")
    filen = "starting_points.txt"
    filp = os.path.join(path1,filen)
    start = open(filp,'r')
    import ast
    xxx = ast.literal_eval(start.readlines()[it])
    start.close()

    xxx = nm.array(xxx)
    xc = nm.array(x)
    print(xxx)
    print(xc)
    xa = xc*xxx
    print(xa)
    # //

    # //
    # Here the procedure file for MARC is created with the "material3d" function in the app.py file.
    # In case a Mooney-Rivlin two parameter model, the correct procedure file will be created.
    # This example is for the Mooney-Rivlin three parameter model
    time.sleep(1)
    if len(x)==3:
        material3d(xa)     #creates the procedure file for Marc
    elif len(x)==2:
        material2d(xa)
    # //

    #//
    # With a pipeline this big and different files and function being imported, a fail safe is added 
    # to sleep the code for a second to ensure it doesn't maybe skip a line or run to fast before 
    # the data from the previous function has been obtained.
    # The subprocess function is used to run the MARC analysis in the background, therefore pausing 
    # the code until MARC is closed.
    time.sleep(1)
    filem = "mat.proc"
    p = subprocess.Popen([mentat_version,filem], bufsize=2048)  #start MARC and load the procedure file 
    # which will open the correct FE model and change the material properties, start MARC solver and 
    # to save the post file for the current DOT increment, close MARC and continue with the code below
    p.wait()
    time.sleep(5)   # the time delay allows Marc to create all output files, even after the program has closed, 
    # it ensures that the code won't crash by loading all the output files too early
    # //
    
    # //
    # Define if MARC analysis converged by loading the NUMERICAL model's status file and reading which
    # EXIT CODE was obtained.
    # - 3004: Converged analysis
    # - Anything other than 3004, no convergence
    sts = filepaths("fem_sts")
    conver = pd.read_csv(sts, header=None, skipinitialspace=True, sep=' ', skiprows=4, keep_default_na=False)   
    time.sleep(1)

    if len(conver)>1:
        c = int(conver.iloc[-3,6])
    else:
        conver = pd.read_csv(sts, header=None, skipinitialspace=True, sep=' ', keep_default_na=False)
        c = int(conver.iloc[-3,6])

    # the boolean constraints are specified
    if c == 3004:
        g[0] = -1       # 3004 says the NUMERICAL model converged and the constraint is satisfied
    else:
        g[0] = 1        
        if len(x)==3:
            xv = [xa[0], xa[1], xa[2]]
        elif len(x)==2:
            xv = [xa[0], xa[1]]
        
        # Here the violated design variables are saved into a file.
        time.sleep(1)
        from functions import violated_constr
        violated_constr(xv,c)     # This function saves the parameters which caused non-convergence 
        # and also what the exit number was
    print(g[0])
    # //

    # //
    # List of all the output and input files needed to create the full data fields.
    # Here the RBF function is called.
    from RBF import RBF_int
    fname1,fname2,fname3,fname4,fnamef,fnamee = filepaths("fem_out","pointfem","exp_out","pointexp","fem_dat","exp_dat")

    ce,de,fci,fdi,fvr,nne,nnf,dm = RBF_int(fname1,fname2,fname3,fname4,fnamef,fnamee,g[0]) # calling the 
    # RBF_int function from the RBF file,
    # - it fills the data of the NUMERICAL model to match that of the EXPERIMENTAL model
    time.sleep(1)
    from RMS import RMS_disp
    rms,rmsd = RMS_disp(de,fvr) # The RMS function is loaded from the RMS file to determine the RMS value for convergence
    time.sleep(1)
    print ("RMS value = ",rms)

    obj.value = rms     #set the objection function of DOT to the RMS value 

#--------------------------------------------------------------------

nDvar = 3       # Tell DOT how many variables need to be solved
if nDvar <=3:
    nCons = 1       #  Tell DOT how many contraint equations are present
else:
    nCons = 2       #  Tell DOT how many contraint equations are present
NEWITR = 0      # Set the current iteration to zero
JWRITE = 21     # File number to write iteration history to.
IWRITE = 20     # File number for printed output


FDCH = 0.001
FDCHM = 0.0001

# // Obtain 10 different design points (starting points) from LHC function
from LHC import LHC
s = 10   #number of different starting points
sp,cl,cu = LHC(nDvar,s)
# //

# ----------------------------------------------------------------------------------------------------------------
# In case more than one algorithm needs to be tested, but the same starting points. Unindent the following lines and
# manually add the starting points from the "starting_points.txt" file. This is to be done before the second algorithm
# is started. During the first, leave the file as is and for the 2nd, unindent the lines below and follow the steps 
# set out below.

# NB !!!! - Remember to indent lines 173-175
#         - Change the algorithm type in line 273
#         - In the "graphs.py" file unindent the lines explained there, only just before the second algorithm's 
#           optimisation is started
#         - Change the lines as explained in "File_paths.py"

# sp = nm.empty([10,3])
# sp[0,:] = [0.23985539715563783, 0.1086356159032091, 0.06752244405180519]
# sp[1,:] = [0.22346941962260952, 0.08517834975399942, 0.04888277372798288]
# sp[2,:] = [0.23680911902278512, 0.1008993817417827, 0.05522682051951559]
# sp[3,:] = [0.2667989705559038, 0.09145828980269591, 0.047717784644901534]
# sp[4,:] = [0.2780371475302363, 0.08052791834324839, 0.06226274627149697]
# sp[5,:] = [0.21648774294463674, 0.10232477181284882, 0.0613617606707654]
# sp[6,:] = [0.2846377298682449, 0.11539898192753328, 0.05897444807215557]
# sp[7,:] = [0.2564777729233578, 0.11022702843829285, 0.05185050155523324]
# sp[8,:] = [0.3009270768118156, 0.08863017446473127, 0.0659901485250131]
# sp[9,:] = [0.3083142113950629, 0.09628959633853354, 0.05326192489971376]
# cl = [0.208454096, 0.078039848, 0.046000552]
# cu = [0.312681144, 0.117059772, 0.069000828]
# ------------------------------------------------------------------------------------------------------------------

# // the next few lines of code creates a new folder containing the Results and set the path for the file as well
path1 = filepaths("path1")
access_rights = 0o777
if not os.path.exists(path1):
    os.mkdir(path1,access_rights)
# //

# // Next a new file is created in the new Results folder containing the list of starting points from LHC as well
# as the limits for the starting points.
# The correct path to the folder is created and given for the file
filename = "starting_points.txt"  
filepath = os.path.join(path1,filename) 

st = open(filepath,"w")       # the starting points file is written 
for sl in range(len(sp)):
    st.write("%s\n"%list(sp[sl]))
lines_w = ["%s\n"%cl, "%s"%cu]
st.writelines(lines_w)       # the writelines command writes all the strings in the list, lines_w, in its own line
st.close()
# //

for ii in range(len(sp)):  # This for loop will allow the program to run and optimise for each starting point
    x  = sp[ii]/sp[ii]   # First starting point in list
    xl = cl/sp[ii]     # Lower bound (lower side constraint) same as LHC
    xu = cu/sp[ii]     # Upper bound (upper side constraint) same as LHC 
    
    # // Create a file which will save the current starting points optimisation number (1 or 2 or 3 ..), for when
    # DOT needs to write out the results it will use the correct file
    iters = open("iterations.txt","w")
    l = "%s"%ii
    iters.write(l)
    iters.close()
    #//

    # // Write the starting point to a results file, this same file will be used after the optimisation to write
    #  the end results.
    filename1 = "results%s.py"%ii
    filepath1 = os.path.join(path1,filename1)     # create the file in the Results folder, and set the path to
    # this file and folder.  
    # //

    # //
    # Create the violated constraints file
    filename2 = "violated%s.txt"%ii
    filepath2 = os.path.join(path1,filename2)    # create the file in the Results folder, and set the path to 
    # this file and folder.
    # Write comments in the "violated constraints" file
    v = open(filepath2,"w") 
    v.write("Material properties resulted in no convergence: \n")
    v.close()
    # //

    # //
    # Write comments in the "results" file.
    # - the starting design point is written to this file. After optimisation the final design point will written 
    #   to this file as well.
    o = open(filepath1,"w")
    o.write("import matplotlib.pyplot as plt")
    o.write("#starting values \n")
    st = "xs = %s \n"%list(sp[ii])          # convert the starting point ndarray to a list
    o.write(st)
    o.close()       # close the file for now so that it can be opened later to append the rest of the data in 
                    # another piece of code
    # // 

    aDot = dot.dot()

    aDot.nPrint    = 5      #Print only the final objection function value, the contraint final values and the final design variables value
    aDot.nMethod   = 3      #Use MMFD as optimisation = 1, SQP = 3, SLP = 2
    aDot.nmRPRM[8] = FDCH
    aDot.nmRPRM[9] = FDCHM
    aDot.nmIPRM[4] = IWRITE  
    aDot.nmIPRM[12] = JWRITE 
    aDot.evaluate  = myEvaluate
    print(aDot.dotcall( x, xl, xu, nCons ))     #call DOT

# To draw some graphs for results, import the following function
from graphs import post_process
# 2 or 3 depending on the dimension
best = post_process(3)
time.sleep(1)
