# -----------------------------------------------------------------------------------------------------------------
# NB: This file was written for a specific project to use linear and RBF interpolation schemes. This file can be rewritten
# as desired by the user, but the rest of the pipeline will have to be adjusted accordingly.

# This file was coded by assuming that the FE analysis used remeshing. Therefore the input file is used to obtain the
# nodes for tracking. The code is written by assuming the FE models were constructed in a specific sequence and
# the instructions to the whole pipeline needs to be read first.

# Francè Bresler
# 12 February 2020
# -----------------------------------------------------------------------------------------------------------------

import numpy as np 
import pandas as pd 
import os
from py_post import *
from py_mentat import *
import scipy.interpolate as sp
import os
import csv
import time
from functions import read_dat_node_number, read_dat_element_number, get_node_position, rbf_interp
from code_settings import pipeline_files
pf = pipeline_files() 

# // In this file the RBF_int function interpolates the nodal data from the "NUM" model to the same nodal locations of the
# "EXP"  model to obtain the same number nodes for both data sets. This is done by first obtaining the locations of each
# node at each increment for the y-displacement of the indentor in both the "NUM" and "EXP" data. Since
# the "NUM" data needs to fit the "EXP" data, the "EXP" data's node coordinates are used. Also both "NUM" and "EXP" data are
# obtained through Marc which is an iterative solver, therefore at each solving increment of the indentor's y-displacement, there
# are more than one sub-increment with the last sub-increment at that solver increment is the sub-increment which contains the
# converged data for the level,thus only the last sub-incremental data are used. 
# After the data has been organised into the indentor increments, for the specific number of nodes presents in the "NUM"
# and "Exp" data respectively, the interpolation can start. First the "NUM" data is linearly interpolated to the missing increments 
# for the "EXP" data's indentor displacements. There after Radial basis functions (RBF) are used to interpolate the current
# level's "NUM" data to the nodal locations of the "EXP" data.
#   -   fname1 = NUM output file (.t16 file)
#   -   fname2 = NUM indentor level/increment output file
#   -   fname3 = EXP output file (.t16 file)
#   -   fname4 = EXP indentor level/increment output file
#   -   fnamef = NUM input file (.dat file)
#   -   fnamee = EXP input file (.dat file)
#   -   G = convergence code, -1 for converged or 1 for no convergence


def RBF_int(fname1, fname2, fname3, fname4, fnamef, fnamee, G):
    # Since remeshing is used with tet10 elements, tracking of nodes are not possible. Therefore the membrane element nodes are tracked
    # Also the output file gave error with the py_mentat API and the nodes present in he sets in the output file was not always obtainable
    # but it is present in the input file and that is why the input files are loaded.
    # Import NUM input file
    datf_nm  = pd.read_csv(fnamef,header=None, names=list(range(27)), skipinitialspace=True, sep=' ', keep_default_na=False)
    # Import EXP input file
    datf_exp = pd.read_csv(fnamee,header=None, names=list(range(27)), skipinitialspace=True, sep=' ', keep_default_na=False)
    femd = post_open(fname1) # open the post file for the NUM data
    expd = post_open(fname3) # open the post file for the EXP data
    
    # Saved the indentors displacement increments during the simulations and the next few lines will import those y-displacement locations files
    location_fem_points_file = pd.read_csv(fname2,header =None, names =['X','Y'], skipinitialspace=True, skiprows=[0,1,2,3,4,5,6,7,8], sep=' ', keep_default_na=False)
    location_points_file = pd.read_csv(fname4,header =None, names =['X','Y'], skipinitialspace=True, skiprows=[0,1,2,3,4,5,6,7,8], sep=' ', keep_default_na=False)
    femd.moveto(1) # move to the first increment to obtain all relevant data and not just the model data.
    expd.moveto(1)

    fnns = femd.node_scalars()      # obtain the number nodal scalar values inside the NUM data post file
    fnes = femd.element_scalars()   # obtain the number of element scalars inside the NUM data post file
    fnince = femd.increments()      # obtain the number of increments in the post file of the FEM data
    sets = femd.sets()              # number of sets present in NUM file
    enns = expd.node_scalars()      # obtain the number nodal scalar values inside the EXP data post file
    enes = expd.element_scalars()   # obtain the number of element scalars inside the EXP data post file
    enince = expd.increments()      # obtain the number of increments in the post file of the EXP data
    sets_e = expd.sets()            # number of sets present in EXP file

    # save the locations of the indentor to a list for the EXP data since these will be the locations to interpolate to
    y_disp_load = list(abs(location_points_file.iloc[:,1]))  

    # save the current y-displacement locations of the NUM indentor 
    y_disp_load_fem = list(abs(location_fem_points_file.iloc[:,1]))

    #   It is neccessary to know if the simulation converged, if it did G = -1, and the y-disp for the NUM data will
    #   cover the whole displacemnet range of the indentor, with the last value being the same as the EXP data's last y-disp value.
    #   If it did not converge the NUM data will not cover the whole displacement range and the therefore won't be able to interpolate 
    #   to the maximum displacement, since there is no data. In this case the EXP data will only be used up until the indentor
    #   level the NUM data did converge at. The data is not accurate for the non-convergence instance, which is okay, since it will
    #   produce a bigger objective value therefore adding to the fact that this point will not cause the optimisation to converge
    #   and the optimisation algorithm knows this point violated the constraints and it will search for a new solution.
    if G == -1:
        y_disp_load = y_disp_load
        print("Correct Exit Code")
    else:
        print("Nooo")
        # Determine if the last level in the NUM data did converge at is a possible level in the EXP data
        if len(y_disp_load)>len(y_disp_load_fem):
            last = np.where(np.array(y_disp_load) == y_disp_load_fem[-1])[0]

            if len(last)>0: # if the level is somewhere in the EXP data, the data up until that point will be used and the rest neglected.
                for popn in range((len(y_disp_load)-last[-1])-1):
                    y_disp_load.pop(-1)
            elif len(last)==0:  # if the last level is not present in the EXP data the following will be done:
                lp = len(y_disp_load_fem)   # the number of NUM y-disp level is determined
                for popn in range(len(y_disp_load)-lp): # the EXP data for the same number of level are used and the rest neglected
                    y_disp_load.pop(-1)
                if (y_disp_load_fem[-1] > y_disp_load[-2]): # from the levels left, if the last level in the NUM data is bigger than the second last
                    y_disp_load[-1] = y_disp_load_fem[-1]   # the last level in the EXP data is set equal to the NUM data's last level
                else:
                    # if the last point in the NUM data is still smaller than the last point in the EXP data
                    # and also smaller than the second last point in the EXP data, then the second last point in the Exp
                    # data is set to the NUM last data point and the last point in the EXP data neglected
                    y_disp_load[-2] = y_disp_load_fem[-1]
                    y_disp_load.pop(-1)
        else:
            lp = len(y_disp_load_fem)-len(y_disp_load)
            for popn in range(lp):
                y_disp_load_fem.pop(-1)
            y_disp_load_fem[-1] = y_disp_load[-1]

    # // Here the it is determined if the simulastion is 2D or 3D and it saves the node_scalar location in the output to allow for further referencing
    stress_strain = []
    stress_strainf = []
    columns = []
    #return a list of the dimensions by ID start from 0
    scalars_of_femd = [femd.node_scalar_label(index) for index in range(0,fnns)]
    dimension = list(np.arange(0,len([i for i, j in zip(scalars_of_femd, ['Displacement X', 'Displacement Y', 'Displacement Z']) if i == j])))

    # NUM data
    # Since the input file is needed to obtain which nodes and elements to track in the output file, we need to search where the information
    # is in the input file
    node_tracking_surface = pf.nm_track_surface_nodes #name of the membrane surface nodes used for tracking
    elemt_traking_surface = pf.nm_track_surface       #name of the membrane surface elements used for tracking

    nn_num,nnseqf = read_dat_node_number(fnamef,femd,node_tracking_surface,'attach')    #node information
    nsq = read_dat_element_number(fnamef,femd,elemt_traking_surface)   #element information

    # // The data arrays are created to split and store the data accordingly. 3D arrays are created:
    #   - the first input is the number of iterations in the simulation (not iterative results)
    #   - the second input is the number of nodes as the rows
    #   - the third input is the dimension as columns, (or the stress & strain components)
    #   - fem_str contains the node number in column 0 and the strain & stress components in the "stress_strainf" list 
    #   from column 1 onwards. 

    #'get_node_position' function extracts all information regarding the nodal positions and displacements for all time steps
    fem_coor,fem_disp = get_node_position(femd,nnseqf,range(len(y_disp_load_fem)),len(dimension))

    # EXP data
    # // The next section is the same as the NUM data section but only for the EXP data files
    node_tracking_surface = pf.exp_track_surface_nodes  #name of the membrane surface nodes used for tracking
    elemt_traking_surface = pf.exp_track_surface        #name of the membrane surface elements used for tracking

    nn_exp,nnseq = read_dat_node_number(fnamee,expd,node_tracking_surface,'attach')
    nsq_e = read_dat_element_number(fnamee,expd,elemt_traking_surface)   #element information
    exp_coor,exp_disp = get_node_position(expd,nnseq,range(len(y_disp_load_fem)),len(dimension))
            
    femd.close()
    expd.close()

    #// Interpolate each node linearly in the NUM data to match the Experimental data's indentor levels
    # This loop determines which indenter level is the same for the EXP and NUM data, this data is 
    # then stored as is, the other levels in the NUM data that does not match the EXP data is zero 
    # for the moment and data needs to be interpolated to fill them
    fem_coor_int = np.empty([len(y_disp_load),len(nn_num),len(dimension)])
    fem_disp_int = np.empty([len(y_disp_load),len(nn_num),len(dimension)])

    ind = [y_disp_load_fem[i] == y_disp_load[i] for i in range(len(y_disp_load))]
    ind_nt = [not(item) for item in ind]
    fem_coor_int[ind,:,:] = fem_coor[ind,:,:]
    fem_disp_int[ind,:,:] = fem_disp[ind,:,:]
    test_fem = np.where(ind_nt)

    # // At the levels which were zero the data are now linearly interpolated. The closest two levels from the NUM data are used as the 
    # bounds and the EXP level is the value in between the bound to interpolate to. All the coordinates, displacements, stress and strain
    # data are interpolated to the level for each NUM node
    for i in range(fem_coor.shape[1]):
        
        fem_coor_int[:,i,0] = sp.interp1d(y_disp_load, fem_coor[:,i,0])(y_disp_load_fem)
        fem_coor_int[:,i,1] = sp.interp1d(y_disp_load, fem_coor[:,i,1])(y_disp_load_fem)
        fem_coor_int[:,i,2] = sp.interp1d(y_disp_load, fem_coor[:,i,2])(y_disp_load_fem)
        
        fem_disp_int[:,i,0] = sp.interp1d(y_disp_load, fem_disp[:,i,0])(y_disp_load_fem)
        fem_disp_int[:,i,1] = sp.interp1d(y_disp_load, fem_disp[:,i,1])(y_disp_load_fem)
        fem_disp_int[:,i,2] = sp.interp1d(y_disp_load, fem_disp[:,i,2])(y_disp_load_fem)
        
    # The rbf_interp function below interpolates data passed to it so that the experimental and fem data sets are interpolated to have the same
    # number of entries, this allows for direct comparison 
    
    interpolated_data = rbf_interp(fem_coor_int,exp_coor,fem_disp_int,pf.interp_direc)

    print(interpolated_data.shape)
    print("RBF Phase Complete")
 
    # // Below all the information needed is written out:
    # EXP coordinates (x,y and/or z), EXP displacements (x,y and/or z), NUM linear interpolated coordinates (2D/3D),
    #  NUM linear interpolated displacements (2D/3D), NUM RBF interpolated data (displacements (2D/3D), number of EXP nodes,
    # number of NUM nodes, dimension of data (2D,3D)
    # - abbreviations of order above: ce,de,fci,fdi,fvr,nn_exp,nnf,dm
    
    return(exp_coor,exp_disp,fem_coor_int,fem_disp_int,interpolated_data,len(nn_exp),len(nn_num),len(dimension))