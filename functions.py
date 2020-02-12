# ---------------------------------------------------------------------------------------------------------------

#  This folder contains all the functions which either writes files/folders or saves data to specific folders
# Each function is explained accoringly

# Francè Bresler
# 12 February 2020

# ---------------------------------------------------------------------------------------------------------------

import numpy as np
import os
import shutil
from File_paths import filepaths

# // The append_val function loads the result file created for the specific optimisations starting point and 
# saves the final data after optimisation in the same file
# - It takes the inputs, X = as a list with the 2/3 optimal coefficients, n_iter = list of the number of iterations,
#   n_obj = list of the objective function values at each iteration
def append_val(X,n_iter,n_obj):
    file = open("iterations.txt","r")   # loads the iterations file to determine which starting points data this is
    it = int(file.readline())
    file.close

    filen = "results%s.py"%it   # the file where the final results will be written to, the path1 line below changes 
    # for each NUMERICAL model in the "File_paths" function
    
    # // The path1 line below needs to change for each simulation or depending where the files needs to be saved to
    path1 = filepaths("path1")
    filepath = os.path.join(path1,filen)
    if not os.path.exists(path1):
        os.makedirs(path1)

    with open(filepath, "a+") as c:
        c.write("#final values \n")
        st = "xf = %s \n"%X 
        si = "iterations = %s \n"%n_iter        # the results file is a .py file which only stores the begin and end point,
        sb = "objective = %s \n"%n_obj          # the number of iterations, and the objective function at each iteration, 
        c.write(st)                             # it also plot the results to show the convergence
        c.write("import matplotlib.pyplot as plt \n")
        c.write(si)
        c.write(sb)
        c.write("plt.plot(iterations,objective,'r-') \n")
        c.write("plt.xlabel('iteration') \n")
        c.write("plt.ylabel('Objective') \n")
        c.write("plt.show() \n")
        c.close()
# //

#// This function opens the violated text file, this function is called everytime a simulation did not converge 
# and stores the points.
# This function takes the input, X = a list of the coefficients, con = the exit number written out by Marc
def violated_constr(X,con):
    file = open("iterations.txt","r")   # Calls the iterations file to know during which strating points optimisation
                                        # the constraints were violated
    it = int(file.readline())
    file.close

    filen = "violated%s.txt"%it

    #// The path1 line needs to change for each simulation, depending where the file is stored
    path1 = filepaths("path1")
    filepath = os.path.join(path1,filen)
    if not os.path.exists(path1):
        os.makedirs(path1)

    with open(filepath, "a+") as v:
        vs = "x = %s \n"%X                  # the function only open the violated file to store the coefficient
                                            # which caused non-convergence
        vc = "Exit number: %s \n"%con       # It also saves the exit number given by Marc
        v.write("\n")
        v.write("Material properties \n")
        v.write(vs)
        v.write(vc)
        v.close()
# //

# // The material function, creates a new procedure file with the updated varaibles.
# the function only takes the optimal point as a list X.
# - material3d: Mooney-Rivlin three parameter model (C10, C01, C20)
# - material2d: Mooney-Rivlin two parameter model (C10, C01)
def material3d(x):
    matp = filepaths("mat_app_path")
    with open('mat.proc', "w") as c:
        c.write(matp)
        c.write('\n')
        c.write('*edit_mater rubber\n')
        c.write('\n')
        m4 = '*mater_option structural:type:mooney\n'
        m5 = '*mater_option structural:mooney_model:five_term\n'
        m1 = '*mater_param structural:mooney_c10 %s\n'%x[0]
        m2 = '*mater_param structural:mooney_c01 %s\n'%x[1]
        m3 = '*mater_param structural:mooney_c20 %s\n'%x[2]
        c.write(m4)
        c.write(m5)
        c.write(m1)
        c.write(m2)
        c.write(m3)
        c.write('*edit_job job1\n')
        c.write('*update_job\n')
        c.write('*save_model\n')
        c.write('*submit_job 1 *monitor_job\n')
        c.write('@popdown(job_properties_pm) @main(results) @popup(modelplot_pm) *post_open_default\n')
        c.write('*save_model\n')
        c.write('*history_collect 0 999999999 1\n')
        c.write('*prog_option history_plot:data_carrier_type_y:cbody\n')
        c.write('*set_history_data_carrier_cbody_y\n')
        c.write('load\n')
        c.write('*set_history_cbody_variable_y\n')
        c.write('Pos Y\n')
        c.write('*history_add_curve\n')
        c.write('*history_fit\n')
        c.write('*history_write pointfem yes\n')
        c.write('*post_close\n')
        c.write('*save_model *quit yes')
        c.close()

def material2d(x):
    matp = filepaths("mat_app_path")
    with open('mat.proc', "w") as c:
        c.write(matp)
        c.write('\n')
        c.write('*edit_mater rubber\n')
        c.write('\n')
        m4 = '*mater_option structural:type:mooney\n'
        m1 = '*mater_param structural:mooney_c10 %s\n'%x[0]
        m2 = '*mater_param structural:mooney_c01 %s\n'%x[1]
        c.write(m1)
        c.write(m2)
        c.write('*submit_job 1 *monitor_job\n')
        c.write('@popdown(job_properties_pm) @main(results) @popup(modelplot_pm) *post_open_default\n')
        c.write('*history_collect 0 999999999 1\n')
        c.write('*prog_option history_plot:data_carrier_type_y:cbody\n')
        c.write('*set_history_data_carrier_cbody_y\n')
        c.write('load\n')
        c.write('*set_history_cbody_variable_y\n')
        c.write('Pos Y\n')
        c.write('*history_add_curve\n')
        c.write('*history_fit\n')
        c.write('*history_write points1 yes\n')
        c.write('*post_close\n')
        c.write('*save_model *quit yes')
        c.close()
# // 

# NB : the "expdata" and "fem_orig_data" functions will be better understood after going through the "RBF" file.
# // This function takes the interpolated data from the RBF function and stores the final data obtained by the optimisation.
#   The function takes the inputs:
#   ce = exp data coordinates
#   de = exp data displacement

#   fvr =  Numerical interpolated data: displacements and Von Mises total strain & stress data or other data captured by user.
#   nne = number of nodes in the experimental data
#   *args = if nothing is specified here the default is to read the iteration number from the iterations.txt file.
#           - Only one argument further needed, the file to save the results to and the iteration number will
#             automatically be equal to 0.
def expdata(ce,de,fvr,nne,*args):
    n_inc = ce.shape[0]
    cc = ce.shape[-1]
    dc = de.shape[-1]
    col = cc+dc
    if len(args)==0:
        file = open("iterations.txt","r")
        it = int(file.readline())
        file.close
        
        filen = "exp_data%s.npy"%it
        filen1 = "RBF_data%s.npy"%it
        path1 = filepaths("path1")

    elif len(args)>0:
        it = 0
        filen = "exp_data%s.npy"%it         # "ce", "de" and "stot" is saved to this file
        filen1 = "RBF_data%s.npy"%it        # "fvr" is saved to this file
        pat = filepaths("path1")
        fl = args[0]
        path1 = os.path.join(pat,fl)

    filepath = os.path.join(path1,filen)
    filepath1 = os.path.join(path1,filen1)
    if not os.path.exists(path1):
        os.makedirs(path1)
    
    expdata = np.empty([n_inc,nne,col])
    expdata[:,:,:cc] = ce
    expdata[:,:,cc:col] = de

    np.save(filepath,expdata)   # save experimental data to this file
    np.save(filepath1,fvr)      # save the numerical interpolated data to this file
    # The data saved separately as python .npy files for easy access later. 
# //

# // This function takes the ORIGINAL NUMERICAL data before interpolation from the RBF function and stores the final data obtained by the optimisation
#    the function takes the inputs:
#   fci = original NUMERICAL data coordinates
#   fdi = original NUMERICAL data displacement
#   nnf = number of nodes in NUMERICAL data
#   nne = number of nodes in the experimental data
#   *args = if nothing is specified here the default is to read the iteration number from the iterations.txt file
#           - Only one argument further needed, the file to save the results to and the iteration number will
#             automatically be equal to 0
def fem_orig_data(fci,fdi,nnf,nne,*args):
    if len(args)==0:
        file = open("iterations.txt","r")
        it = int(file.readline())
        file.close
        filen = "fem_orig_data%s.npy"%it
        path1 = filepaths("path1")

    elif len(args)>0:
        it = 0
        filen = "fem_orig_data%s.npy"%it    # "fci", "fdi" is saved to this file
        pat = filepaths("path1")
        fl = args[0]
        path1 = os.path.join(pat,fl)

    filepath = os.path.join(path1,filen)
    if not os.path.exists(path1):
        os.makedirs(path1)
    
    n_inc = fci.shape[0]
    fcd = fci.shape[-1]
    fdc = fdi.shape[-1]
    col = fcd + fdc
    femdata = np.empty([n_inc,nnf,col])
    femdata[:,:,:fcd] = fci
    femdata[:,:,fcd:col] = fdi

    np.save(filepath,femdata)   # saves the original FEM data as a .npy file
# //

# // This function only stores the final values obtained after optimisation for each starting point.
#   It takes the input X = the optimal point, 
#   *args = if none is given, the file will be saved to the results folder.
#           - Only 1 input needed, the folder where to store the file
def final_points(X,*args):
    filen = "final_points.txt"
    if len(args)==0:
        path1 = filepaths("path1")
    elif len(args)>0:
        pat = filepaths("path1")
        fl = args[0]
        path1 = os.path.join(pat,fl)
          
    filepath = os.path.join(path1,filen)
    if not os.path.exists(path1):
        os.makedirs(path1)

    with open(filepath, "a+") as c:
        st = "%s \n"%X 
        c.write(st)
        c.close()
# //

#// This function creates a new procedure file which changes the:
#    X = material coefficients (Mooney-Rivlin, three parameter (C10, C01, C20))
#    fm = the remesh edge length.
def remesh_proc(x,fm):
    matp = filepaths("mat_app_path")
    with open('remesh.proc', "w") as c:
        c.write(matp)
        c.write('\n')
        c.write('*edit_mater rubber\n')
        c.write('\n')
        m1 = '*mater_param structural:mooney_c10 %s\n'%x[0]
        m2 = '*mater_param structural:mooney_c01 %s\n'%x[1]
        m3 = '*mater_param structural:mooney_c20 %s\n'%x[2]
        c.write(m1)
        c.write(m2)
        c.write(m3)
        c.write('*edit_adapg adapg1\n')
        fine = '*adapg_dens_ctrl_param 1 target_edge_length %s\n'%fm
        fine1 = '*adapg_dens_ctrl_param 2 small_edge_length_cv %s\n'%fm 
        c.write(fine)
        c.write(fine1)
        c.write('*edit_job job1\n')
        c.write('*update_job\n')
        c.write('*save_model\n')
        c.write('*submit_job 1 *monitor_job\n')
        c.write('@popdown(job_properties_pm) @main(results) @popup(modelplot_pm) *post_open_default\n')
        c.write('*save_model\n')
        c.write('*history_collect 0 999999999 1\n')
        c.write('*prog_option history_plot:data_carrier_type_y:cbody\n')
        c.write('*set_history_data_carrier_cbody_y\n')
        c.write('load\n')
        c.write('*set_history_cbody_variable_y\n')
        c.write('Pos Y\n')
        c.write('*history_add_curve\n')
        c.write('*history_fit\n')
        c.write('*history_write pointfem yes\n')
        c.write('*post_close\n')
        c.write('*save_model *quit yes')
        c.close()
# //

# // This function saves for each optimisation procedure corresponding to the starting point from the LHC function,
#   the objective function value from each optimisation iteration and its corresponding iteration number.
#   obj = list of objective function values
#   ite = list of iteration numbers
#   - They are saved as .npy files for eacy access during post processing.
def objectfunc(obj,ite):
    file = open("iterations.txt","r")
    it = int(file.readline())
    file.close
    objf = "objectfunction%s.npy"%it    # "obj" is saved to this file
    iterf = "niter%s.npy"%it            # "ite" is saved to this file
    path1 = filepaths("path1")

    objpat = os.path.join(path1,objf)
    iterpat = os.path.join(path1,iterf)
    if not os.path.exists(path1):
        os.makedirs(path1)

    np.save(objpat,obj)
    np.save(iterpat,ite)
# //
    