# All filepaths needed in the whole pipeline.
#-------------------------------------------------------------------------------------------------------------
# This file contains the necessary filepaths to store and access the necessary files and data saved during the 
# whole pipeline.
# The user just add their own file path to the "InverseMatMarc" folder and leave the rest as is.

# More than one optimisation algorithm  can be used
# Lines 59 -75 explain what needs to change in this file during such a procedure

# Francè Bresler
# 12 February 2020
#-------------------------------------------------------------------------------------------------------------

import numpy as np

# ------------------------------------------------------------------------------------------------------------
# The function works by assigning it to a variable and calling the desired path variable by using the name given.
# More than one path can be accessed as a tuple.
# ex. > from File_paths import filepaths
#     > path1 = filepaths("fem_out")
#     OR
#     > path1, path2 = filepaths("fem_out", "pointfem")

def filepaths(*args):
    from code_settings import pipeline_files
    pf = pipeline_files()
    output = []
    for i in range(len(args)):
        # Marc's .t16 output file for the "Numerical" FE model 
        if (args[i]) == "fem_out":      # Command to call with function to obtain the filepath specified below
            output.append(pf.crnt_fldr("fwd") + '/' + pf.sub_folder + '/'+ pf.nm_t16)
        # The file created from the "Numerical" simulation, which contains the indentation level for each increment
        elif (args[i]) == "pointfem":
            output.append(pf.crnt_fldr("fwd") + '/'+ "pointfem")
        # Marc's .t16 output file for the "Experimental" FE model
        elif (args[i]) == "exp_out":
            output.append(pf.crnt_fldr("fwd") + '/' + pf.sub_folder + '/'+ pf.exp_t16)
        # The file created from the "Experimental" simulation, which contains the indentation level for each increment
        elif (args[i]) == "pointexp":
            output.append(pf.crnt_fldr("fwd") + '/' + pf.sub_folder + '/'+ "cyl_diag/pointexp")
        # The input file for the "Numerical" model
        elif (args[i]) == "fem_dat":
            output.append(pf.crnt_fldr("fwd") + '/' + pf.sub_folder + '/'+ pf.nm_dat)
        # The input file for the "Experimental" model
        elif (args[i]) == "exp_dat":
            output.append(pf.crnt_fldr("fwd") + '/' + pf.sub_folder + '/'+ pf.exp_dat)
        # The status file for the "Numerical" model, to obtain the EXIT CODE
        elif (args[i]) == "fem_sts":
            output.append(pf.crnt_fldr("fwd") + '/' + pf.sub_folder + '/'+ pf.nm_sts)
        # Filepath to store and obtain the procedure file to control the "Numerical" model simulations
        elif (args[i]) == "mat_proc_path":
            output.append(pf.crnt_fldr("fwd") + '/'+ "mat.proc")
        # The name of the procedure file, incase it needs changing and don't want to go through the whole pipeline (OPTIONAL)
        elif (args[i]) == "mat_proc_file":
            output.append("mat.proc")
        # Filtepath and folder where to OBTAIN / STORE the data after DOT completed the optimisation
        elif (args[i]) == "path1":
            output.append(pf.crnt_fldr("fwd") + '/' + pf.sub_folder + '/'+ "ResultsSQP")
        # Marc Mentat version 2019    
        elif (args[i]) == "2019":
            output.append("C:/Program Files/MSC.Software/Marc/2019.0.0/mentat2019/binmentat.bat")
        # Marc Mentat version 2020
        elif (args[i]) == "2020":
            output.append("C:/Program Files/MSC.Software/Marc/2020.0.0/mentat2020/bin/mentat.bat")  

# -------------------------------------------------------------------------------------------------------------------        
        # In case two optimisation algorithms are tested:
        # - During the first algorithm optimisation name the folder as desired; here is was named "ResultsSQP" 
        #   with the SQP part refering to the optimisation algorithm.
        # - Before going to the next algorithm, set path1 to second algorithm's results and set path2 to the 
        #   first algorithm's results.
        # - ex. 1st algorithm run:
        #        > path1 = ..../ResultsSQP
        #       2nd algorithm run:
        #        > path1 = ..../ResultsSLP
        #        > path2 = ..../ResultsSQP

        # NB !!!!!! - Remember to change the main_code file and the graphs file
        elif (args[i]) == "path2":
            output.append(pf.crnt_fldr("fwd") + '/' + pf.sub_folder + '/'+ "ResultsSLP")
        # Filtepath and folder where to OBTAIN / STORE the data after DOT completed the optimisation (OPTIONAL)
        elif (args[i]) == "path3":
            output.append(pf.crnt_fldr("fwd") + '/' + pf.sub_folder)
# -------------------------------------------------------------------------------------------------------------------
        
        # Filepath to access the "Numerical" model, to be used to create the procedure file where MARC will then open this 
        # FE model for analysis.
        elif (args[i]) == "mat_app_path":
            output.append(r'*open_model ' + pf.crnt_fldr("bwd") + '\\' + pf.sub_folder + '\\' + pf.nm_mud)
    if len(output)>1:
        output = tuple(output)
        return(output)
    else:
        return(output[0])




