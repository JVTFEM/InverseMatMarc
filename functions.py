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
        c.write('*model_orientation_iso1\n')
        c.write('identify_contact *redraw\n')
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
        
        exp_data_File = "exp_data%s.npy"%it
        RBF_data_File = "RBF_data%s.npy"%it
        path1 = filepaths("path1")

    elif len(args)>0:
        it = 0
        exp_data_File = "exp_data%s.npy"%it         # "ce", "de" and "stot" is saved to this file
        RBF_data_File = "RBF_data%s.npy"%it         # "fvr" is saved to this file
        pat = filepaths("path1")
        fl = args[0]
        path1 = os.path.join(pat,fl)

    filepath  = os.path.join(path1,exp_data_File)
    filepath1 = os.path.join(path1,RBF_data_File)
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




#this function performs the task of extracting the node numbers as reference numbers for both the .t16 and .dat files
def read_dat_node_number(dat_file,t_16_file,start_row,end_row):
    import pandas as pd
    t_16_file.moveto(1)
    datf_nm  = pd.read_csv(dat_file,header=None, names=list(range(27)), skipinitialspace=True, sep=' ', keep_default_na=False)
    #determine search bounds
    start_row  = datf_nm.index[datf_nm.iloc[:,3] == start_row][0] + 3 #determine row position with tracking surface nodes using membrane's name
    end_row    = datf_nm.index[datf_nm.iloc[:,0] == end_row][0] + 2 #determine row position where the node number for membrane are no longer listed
    n_rows = end_row - start_row    #determine the number of rows where the membrane nodes are listed
    #get node data
    nodef  = pd.read_csv(dat_file,header=None, skiprows= range(start_row), nrows = n_rows, skipinitialspace=True, sep=' ', keep_default_na=False)
    nodef_list_flat = list(filter(('c').__ne__,[item for sublist in nodef.values.tolist() for item in sublist])) #filter out all instances of 'c'
    nodef_list_flat = list(filter(('').__ne__,nodef_list_flat)) #filter out all instances of ''

    node_numbers = [int(item) for item in nodef_list_flat]  #ensure all entires of list are integers
    #get the node sequence from the .t16 file
    node_sequence = []
    for i in range(node_numbers[0],node_numbers[-1]+1):
        node_sequence.append(t_16_file.node_sequence(i))
    
    return node_numbers,node_sequence   




#this function performs the task of extracting the element numbers as reference numbers for both the .t16 and .dat files
def read_dat_element_number(dat_file,t_16_file,start_row):
    import pandas as pd
    t_16_file.moveto(1)
    datf  = pd.read_csv(dat_file,header=None, names=list(range(27)), skipinitialspace=True, sep=' ', keep_default_na=False)
    #determine location of element ID data
    start_row  = datf.index[datf.iloc[:,3] == start_row][0] + 3 #determine row position with tracking surface nodes using membrane's name
    n_rows = 1   #only 1 row is needed get element numbers since they are given as a range
    #get element data
    nodef  = pd.read_csv(dat_file,header=None, skiprows= range(start_row), nrows = n_rows, skipinitialspace=True, sep=' ', keep_default_na=False)
    elem_range = nodef.values.tolist()
    element_sequence = [t_16_file.element_sequence(i) for i in range(elem_range[0][0],elem_range[0][2]+1)]

    return element_sequence   




def get_node_position(t_16_file,node_sequence,time_steps,n_dimensions):
    #'get_node_position' function extracts all information regarding the nodal positions and displacements for all time steps
    import numpy as np

    if type(time_steps) == int:
        n_steps = 1
        time_steps = [time_steps-1]
    else:
        n_steps = len(time_steps)

    node_coor = np.zeros([n_steps,len(node_sequence),n_dimensions])
    node_disp = np.zeros([n_steps,len(node_sequence),n_dimensions])

    n_nodes = len(node_sequence)
    n_nodes_range = range(n_nodes)
    t_16_file.moveto(1)

    for node in range(n_nodes):
        Node_info = str(t_16_file.node(node_sequence[node])).replace("\n",",").replace(" ","").replace("=",",").split(",")
        node_coor[0,node,0] = Node_info[2] # x locations
        node_coor[0,node,1] = Node_info[4] # y locations
        node_coor[0,node,2] = Node_info[6] # z locations

    #make sure all node coordinates are the same for all time steps
    node_coor[:,:,:] = node_coor[0,:,:]


    for time_int in range(n_steps):
        t_16_file.moveto(time_steps[time_int]+1)
        disp_x = [t_16_file.node_scalar(node_sequence[nod],0) for nod in n_nodes_range] #list geneterated for all x nodes at current time step
        disp_y = [t_16_file.node_scalar(node_sequence[nod],1) for nod in n_nodes_range] #list geneterated for all y nodes at current time step
        disp_z = [t_16_file.node_scalar(node_sequence[nod],2) for nod in n_nodes_range] #list geneterated for all z nodes at current time step
        
        node_disp[time_int,:,0] = np.array(disp_x)
        node_disp[time_int,:,1] = np.array(disp_y)
        node_disp[time_int,:,2] = np.array(disp_z)

    node_coor = node_coor + node_disp
    return node_coor,node_disp



def get_element_strain(dat_file,t_16_file,node_tracking_surface,elemt_traking_surface,time_steps):
    import pandas as pd
    import numpy as np
    
    if type(time_steps) == int:
        n_steps = 1
        time_steps = [time_steps-1]
    else:
        n_steps = len(time_steps)    

    t_16_file.moveto(1)
    node_numbers,node_numbers_sequenced = read_dat_node_number(dat_file,t_16_file,node_tracking_surface,'attach')
    element_numbers_sequenced = read_dat_element_number(dat_file,t_16_file,elemt_traking_surface)   #element information

    node_instances = np.zeros([len(node_numbers),2])
    node_instances[:,0] = node_numbers
    for element in element_numbers_sequenced:
        for i in range(0,4):
            nn = t_16_file.element_scalar(element,6)[i].id
            ind = nn - node_numbers[0] - 1
            node_instances[ind,1] += 1


    
    labels_strain = [6,7,8]   #[min_strain,int_strain,max_strain]   
    #labels_stress = [15,16,17]#[min_stress,int_stress,max_stress]
    
    label = labels_strain
    node_ID_prin = np.zeros([len(node_numbers),4])
    node_ID_prin[:,0] = node_numbers
    
    avg_min = np.zeros([len(node_numbers),n_steps])
    avg_int = np.zeros([len(node_numbers),n_steps])
    avg_max = np.zeros([len(node_numbers),n_steps])
    
    for time_int in range(n_steps):
        t_16_file.moveto(time_steps[time_int]+1)
        for prin_ind in range(0,3): 
            for element in element_numbers_sequenced:
                for i in range(0,4):
                    nn = t_16_file.element_scalar(element,label[prin_ind])[i].id
                    ind = nn - node_numbers[0] - 1
                    node_ID_prin[ind,prin_ind+1] += t_16_file.element_scalar(element,label[prin_ind])[i].value

        #This accounts for the occurances of each node regarding reacurance since a reacuring node is added to the
        #the total sum. An average is taken by the following equation average = [x]/len(x)
        node_ID_prin[:,1] = node_ID_prin[:,1]/node_instances[:,1]   
        node_ID_prin[:,2] = node_ID_prin[:,2]/node_instances[:,1]   
        node_ID_prin[:,3] = node_ID_prin[:,3]/node_instances[:,1]   

        avg_min[:,time_int] = node_ID_prin[:,1]
        avg_int[:,time_int] = node_ID_prin[:,2]
        avg_max[:,time_int] = node_ID_prin[:,3]
    
    return node_numbers, avg_min, avg_int,avg_max



def get_element_stress(dat_file,t_16_file,node_tracking_surface,elemt_traking_surface,time_steps):
    import pandas as pd
    import numpy as np
    
    if type(time_steps) == int:
        n_steps = 1
        time_steps = [time_steps-1]
    else:
        n_steps = len(time_steps)

    t_16_file.moveto(1)
    node_numbers,node_numbers_sequenced = read_dat_node_number(dat_file,t_16_file,node_tracking_surface,'attach')
    element_numbers_sequenced = read_dat_element_number(dat_file,t_16_file,elemt_traking_surface)   #element information

    node_instances = np.zeros([len(node_numbers),2])
    node_instances[:,0] = node_numbers
    for element in element_numbers_sequenced:
        for i in range(0,4):
            nn = t_16_file.element_scalar(element,6)[i].id
            ind = nn - node_numbers[0] - 1
            node_instances[ind,1] += 1


    
    #labels_strain = [6,7,8]   #[min_strain,int_strain,max_strain]   
    labels_stress = [15,16,17]#[min_stress,int_stress,max_stress]
    
    label = labels_stress
    node_ID_prin = np.zeros([len(node_numbers),4])
    node_ID_prin[:,0] = node_numbers
    
    avg_min = np.zeros([len(node_numbers),n_steps])
    avg_int = np.zeros([len(node_numbers),n_steps])
    avg_max = np.zeros([len(node_numbers),n_steps])
    
    for time_int in range(n_steps):
        t_16_file.moveto(time_steps[time_int]+1)
        for prin_ind in range(0,3): 
            for element in element_numbers_sequenced:
                for i in range(0,4):
                    nn = t_16_file.element_scalar(element,label[prin_ind])[i].id
                    ind = nn - node_numbers[0] - 1
                    node_ID_prin[ind,prin_ind+1] += t_16_file.element_scalar(element,label[prin_ind])[i].value

        #This accounts for the occurances of each node regarding reacurance since a reacuring node is added to the
        #the total sum. An average is taken by the following equation average = [x]/len(x)
        node_ID_prin[:,1] = node_ID_prin[:,1]/node_instances[:,1]   
        node_ID_prin[:,2] = node_ID_prin[:,2]/node_instances[:,1]   
        node_ID_prin[:,3] = node_ID_prin[:,3]/node_instances[:,1]   

        avg_min[:,time_int] = node_ID_prin[:,1]
        avg_int[:,time_int] = node_ID_prin[:,2]
        avg_max[:,time_int] = node_ID_prin[:,3]
    # calculate the average von mises stress for all nodes
    vm_stress = np.sqrt(((avg_max - avg_min)**2 + (avg_max - avg_int)**2 + (avg_min - avg_int)**2)/2)
    return node_numbers, avg_min, avg_int, avg_max, vm_stress



def rbf_interp(fem_coords,exp_coords,data,interp_direc):
    import scipy.interpolate as sp
    from code_settings import pipeline_files
    pf = pipeline_files()
    interp_method = pf.interp_method 
    
    if (interp_direc == 0):
        ref_coords = fem_coords 
        int_coords = exp_coords
        interpolated_data = np.zeros(int_coords.shape)

    else:
        ref_coords = exp_coords 
        int_coords = fem_coords
        interpolated_data = np.zeros(int_coords.shape)
        

    for i in range(ref_coords.shape[0]):

        rbfi_x = sp.Rbf(ref_coords[i,:,0],ref_coords[i,:,1],ref_coords[i,:,2],data[i,:,0],method = interp_method)
        rbfi_y = sp.Rbf(ref_coords[i,:,0],ref_coords[i,:,1],ref_coords[i,:,2],data[i,:,1],method = interp_method)
        rbfi_z = sp.Rbf(ref_coords[i,:,0],ref_coords[i,:,1],ref_coords[i,:,2],data[i,:,2],method = interp_method)

        interpolated_data[i,:,0] = rbfi_x(int_coords[i,:,0],int_coords[i,:,1],int_coords[i,:,2])
        interpolated_data[i,:,1] = rbfi_y(int_coords[i,:,0],int_coords[i,:,1],int_coords[i,:,2])
        interpolated_data[i,:,2] = rbfi_z(int_coords[i,:,0],int_coords[i,:,1],int_coords[i,:,2])
        
    return interpolated_data 