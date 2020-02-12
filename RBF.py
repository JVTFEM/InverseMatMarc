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
import scipy.interpolate as si
import os
import csv
import time
 

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
    datf = pd.read_csv(fnamef,header=None, names=list(range(27)), skipinitialspace=True, sep=' ', keep_default_na=False)
    # Import EXP input file
    datf2 = pd.read_csv(fnamee,header=None, names=list(range(27)), skipinitialspace=True, sep=' ', keep_default_na=False)
    femd = post_open(fname1) # open the post file for the NUM data
    expd = post_open(fname3) # open the post file for the EXP data

    # Saved the indentors displacement increments during the simulations and the next few lines will import those y-displacement locations files
    location_fem_points_file = pd.read_csv(fname2,header =None, names =['X','Y'], skipinitialspace=True, skiprows=[0,1,2,3,4,5,6,7,8], sep=' ', keep_default_na=False)
    location_points_file = pd.read_csv(fname4,header =None, names =['X','Y'], skipinitialspace=True, skiprows=[0,1,2,3,4,5,6,7,8], sep=' ', keep_default_na=False)
    expd.moveto(1) # move to the first increment to obtain all relevant data and not just the model data.
    femd.moveto(1)

    fnns = femd.node_scalars()      # obtain the number nodal scalar values inside the NUM data post file
    fnes = femd.element_scalars()   # obtain the number of element scalars inside the NUM data post file
    fnince = femd.increments()       # obtain the number of increments in the post file of the FEM data
    sets = femd.sets()              # number of sets present in NUM file
    enns = expd.node_scalars()      # obtain the number nodal scalar values inside the EXP data post file
    enes = expd.element_scalars()   # obtain the number of element scalars inside the EXP data post file
    enince = expd.increments()       # obtain the number of increments in the post file of the EXP data
    sets_e = expd.sets()              # number of sets present in EXP file

    # save the locations of the indentor to a list for the EXP data since these will be the locations to interpolate to
    y_disp_load = []
    for i in range(len(location_points_file)):
        ss = location_points_file.iloc[i,1]
        y_disp_load.append(abs(ss))

    # save the current y-displacement locations of the NUM indentor 
    y_disp_load_fem = []
    for j in range(len(location_fem_points_file)):
        s = location_fem_points_file.iloc[j,1]
        y_disp_load_fem.append(abs(s))

    # It is neccessary to know if the simulation converged, if it did G = -1, and the y-disp for the NUM data will
    #   cover the whole displacemnet range of the indentor, with the last value being the same as the EXP data's last y-disp value.
    #   If it did not converge the NUM data will not cover the whole displacement range and the therefore won't be able to interpolate 
    #   to the maximum displacement, since there is no data. In this case the EXP data will only be used up until the indentor
    #   level the NUM data did converge at. The data is not accurate for the non-convergence instance, which is okay, since it will
    #   produce a bigger objective value therefore adding to the fact that this point will not cause the optimisation to converge
    #   and the optimisation algorithm knows this point violated the constraints and it will search for a new solution.
    if G == -1:
        y_disp_load = y_disp_load
        print("hey")
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
    dimension = []
    stress_strain = []
    stress_strainf = []
    columns = []
    for i in range(0, fnns):
        dd = femd.node_scalar_label(i)
        if dd == "Displacement X":
            nlocx_fem = i
            dimension.append(nlocx_fem)    
        if dd == "Displacement Y":
            nlocy_fem = i
            dimension.append(nlocy_fem)
        if dd == "Displacement Z":
            nlocz_fem = i
            dimension.append(nlocz_fem)
    # //

    # NUM data
    # Since the input file is needed to obtain which nodes and elements to track in the output file, we need to search where the information
    # is in the input file
    nd_ind = 'skin_nodes'
    s = datf.index[datf.iloc[:,3] == nd_ind][0] # search in which line the membrane nodes set is in the input file
    nd_st = "attach"
    stop = datf.index[datf.iloc[:,0] == nd_st][0] # search in which line the membrane nodes data stopped
    col_st = "c"
    stc = datf.columns[datf.iloc[s+1,:] == col_st][0] # search until which column the data is given
    sr = range(0,s+2)   # Python uses zero indexing and the input file does not, therefore add one line for that and also one line to where the data starts
    sto = stop-s        # determine for how many lines the data is present
    # Import only that block of data
    nodef = pd.read_csv(fnamef,header=None, names=list(range(stc)), skiprows=sr, nrows = sto, skipinitialspace=True, sep=' ', keep_default_na=False)
    size = nodef.shape
    nnseqf = []
    nn = []

    # // Save each node to track in a list. Again Python uses zero indexing, therefore the node from Python is referenced as one less
    for i in range(1,size[0]):
        for j in range(size[1]):
            val = nodef.iloc[i,j]
            if val != '':
                nod = int(nodef.iloc[i,j])
                seq = femd.node_sequence(nod)
                nn.append(nod)
                if seq >=0:
                    nnseqf.append(seq)
    # //

    # // Same as the nodes, the elements are being searched for and saved
    elem_ind = 'skin'
    se = datf.index[datf.iloc[:,3] == elem_ind][0]
    elem_st = "define"
    stope = datf.index[datf.iloc[:,0] == elem_st][np.where(datf.index[datf.iloc[:,0] == elem_st] > se)[0][0]]
    sr = range(0,se+2)
    sto = stope-se
    elemf = pd.read_csv(fnamef,header=None, names=list(range(4)), skiprows=sr, nrows = sto, skipinitialspace=True, sep=' ', keep_default_na=False)
    nelem = range(int(elemf.iloc[1,0]),int(elemf.iloc[1,2])+1)
    nsq = []
    # the elements are referenced differently in the output file, therefore by using the "element_sequence" command it can be determined
    # which elements we are actually using.
    for i in nelem:
        nsq.append(femd.element_sequence(i))
    #  //

    # // The data arrays are created to split and store the data accordingly. 3D arrays are created:
    #   - the first input is the number of iterations in the simulation (not iterative results)
    #   - the second input is the number of nodes as the rows
    #   - the third input is the dimension as columns, (or the stress & strain components)
    #   - fem_str contains the node number in column 0 and the strain & stress components in the "stress_strainf" list 
    #   from column 1 onwards. 
    fem_disp = np.zeros([len(y_disp_load_fem),len(nn),len(dimension)])
    fem_coor = np.zeros([len(y_disp_load_fem),len(nn),len(dimension)])
    # print(nn)
    # print(nsq)

    # // the next two for loops collect the node displacement and coordinate data for each increment
    # The ".node()" function gives the nodes original coordinates and it is being stored at the first increment of the coordinate array
    for nod in range(len(nnseqf)):
        cc = str(femd.node(nnseqf[nod])).replace("\n",",").replace(" ","").replace("=",",").split(",")
        for p in range(len(cc)):
            if cc[p] == "x":
                fem_coor[0,nod,0] = cc[p+1]
            if cc[p] == "y":
                fem_coor[0,nod,1] = cc[p+1]
            if cc[p] == "z":
                fem_coor[0,nod,2] = cc[p+1]

    # here ".node_scalar()" takes 2 arguments, 1 = node number, 2 = dimension (x,y,z)
    for m in range(1,len(y_disp_load_fem)+1):
        femd.moveto(m)
        for nod in range(len(nnseqf)):
            for d in dimension:
                disp = femd.node_scalar(nnseqf[nod],d)
                fem_disp[m-1,nod,d] = disp  # the displacement is the total displacement from the first coordinate
                fem_coor[m-1,nod,d] = fem_coor[0,nod,d] + disp  # the current coordinate for the node
    #//

    # EXP data
    # // The next section is the same as the NUM data section but only for the EXP data files
    nd_ind = 'skin_nodes'
    nd_st = "attach"
    s = datf2.index[datf2.iloc[:,3] == nd_ind][0]
    stop = datf2.index[datf2.iloc[:,0] == nd_st][0]
    col_st = "c"
    stc = datf2.columns[datf2.iloc[s+1,:] == col_st][0]
    sr = range(0,s+2)
    sto = stop-s
    nodef_e = pd.read_csv(fnamee,header=None, names=list(range(stc)), skiprows=sr, nrows = sto, skipinitialspace=True, sep=' ', keep_default_na=False)
    size = nodef_e.shape
    nnseq = []
    nne = []

    for i in range(1,size[0]):
        for j in range(size[1]):
            val = nodef_e.iloc[i,j]
            if val != '':
                nod = int(nodef_e.iloc[i,j])
                seq = expd.node_sequence(nod)
                nne.append(nod)
                if seq >=0:
                    nnseq.append(seq)

    elem_ind = 'skin'
    elem_st = "define"
    se = datf2.index[datf2.iloc[:,3] == elem_ind][0]
    stope = datf2.index[datf2.iloc[:,0] == elem_st][np.where(datf2.index[datf2.iloc[:,0] == elem_st] > se)[0][0]]

    sr = range(0,se+2)
    sto = stope-se
    eleme = pd.read_csv(fnamee,header=None, names=list(range(4)), skiprows=sr, nrows = sto, skipinitialspace=True, sep=' ', keep_default_na=False)
    nelem_e = range(int(eleme.iloc[1,0]),int(eleme.iloc[1,2])+1)
    nsq_e = []
    for i in nelem_e:
        nsq_e.append(expd.element_sequence(i))

    exp_disp = np.zeros([len(y_disp_load),len(nne),len(dimension)])
    exp_coor = np.zeros([len(y_disp_load),len(nne),len(dimension)])
    # print(nne)

    for nod in range(len(nnseq)):
        cc = str(expd.node(nnseq[nod])).replace("\n",",").replace(" ","").replace("=",",").split(",")
        for p in range(len(cc)):
            if cc[p] == "x":
                exp_coor[0,nod,0] = cc[p+1]
            if cc[p] == "y":
                exp_coor[0,nod,1] = cc[p+1]
            if cc[p] == "z":
                exp_coor[0,nod,2] = cc[p+1]

    for m in range(1,len(y_disp_load)+1):
        expd.moveto(m)
        for nod in range(len(nnseq)):
            for d in dimension:
                disp = expd.node_scalar(nnseq[nod],d)
                exp_disp[m-1,nod,d] = disp
                exp_coor[m-1,nod,d] = exp_coor[0,nod,d] + disp
    # //

    #// Interpolate each node linearly in the NUM data to match the Experimental data's indentor levels
    # This loop determines which indenter level is the same for the EXP and NUM data, this data is 
    # then stored as is, the other levels in the NUM data that does not match the EXP data is zero 
    # for the moment and data needs to be interpolated to fill them
    fem_coor_int = np.empty([len(y_disp_load),len(nn),len(dimension)])
    fem_disp_int = np.empty([len(y_disp_load),len(nn),len(dimension)])

    test_fem = [5.0]*len(y_disp_load)
    for df in range(len(y_disp_load_fem)):
        for de in range(len(y_disp_load)):
            if y_disp_load_fem[df] == y_disp_load[de]:
                for c in dimension:
                    fem_coor_int[de,:,c] = fem_coor[df,:,c]
                    fem_disp_int[de,:,c] = fem_disp[df,:,c]
                test_fem[de] = y_disp_load[de]
                break
            else:
                pass
    # the test_fem array stored which levels are missing in the NUM data and were set as zero, the level 
    # number is then searched where it is seen as zero
    test_fem = np.where(np.array(test_fem) == 5.0)[0]
    #//

    # // At the levels which were zero the data are now linearly interpolated. The closest two levels from the NUM data are used as the 
    # bounds and the EXP level is the value in between the bound to interpolate to. All the coordinates, displacements, stress and strain
    # data are interpolated to the level for each NUM node
    xl = []
    countn = 0
    nn_fem_int = [len(nn)]*len(y_disp_load)
    for l in range(len(test_fem)):
        di = y_disp_load[test_fem[l]]
        for h in range(len(y_disp_load_fem)):
            if di < y_disp_load_fem[h]:
                for n in range(len(nn)):
                    for dm in dimension:
                        fem_disp_int[test_fem[l],n,dm] = fem_disp[h,n,dm] + (y_disp_load_fem[h] - di)*((fem_disp[h,n,dm] - fem_disp[h-1,n,dm])/(y_disp_load_fem[h] - y_disp_load_fem[h-1]))
                        fem_coor_int[test_fem[l],n,dm]= fem_coor[h,n,dm] + (y_disp_load_fem[h] - di)*((fem_coor[h,n,dm] - fem_coor[h-1,n,dm])/(y_disp_load_fem[h] - y_disp_load_fem[h-1]))
                    if (fem_coor_int[test_fem[l],n,0] == 0.0) and (fem_coor_int[test_fem[l],n,1] == 0.0) and (fem_coor_int[test_fem[l],n,2] == 0.0):
                        countn = countn+1
                nn_fem_int[l] = len(nn) - countn
                break
    #//

    # // Now that there is data at the same indenter level for both the NUM and EXP data, the NUM data points can be interpolated at each
    # level to the same points as the EXP data. This is done by using Radial Basis Functions.
    fem_vlak_rbf = np.empty([len(y_disp_load),len(nne),len(dimension)+2])
    for vlak in range(len(y_disp_load)):
        rbfv_x = si.Rbf(fem_coor_int[vlak,:nn_fem_int[vlak],0],fem_coor_int[vlak,:nn_fem_int[vlak],1],fem_coor_int[vlak,:nn_fem_int[vlak],2],fem_disp_int[vlak,:nn_fem_int[vlak],0], method='cubic')
        vx = rbfv_x(exp_coor[vlak,:len(nne),0],exp_coor[vlak,:len(nne),1],exp_coor[vlak,:len(nne),2])

        rbfv_y = si.Rbf(fem_coor_int[vlak,:nn_fem_int[vlak],0],fem_coor_int[vlak,:nn_fem_int[vlak],1],fem_coor_int[vlak,:nn_fem_int[vlak],2],fem_disp_int[vlak,:nn_fem_int[vlak],1], method='cubic')
        vy = rbfv_y(exp_coor[vlak,:len(nne),0],exp_coor[vlak,:len(nne),1],exp_coor[vlak,:len(nne),2])

        rbfv_z = si.Rbf(fem_coor_int[vlak,:nn_fem_int[vlak],0],fem_coor_int[vlak,:nn_fem_int[vlak],1],fem_coor_int[vlak,:nn_fem_int[vlak],2],fem_disp_int[vlak,:nn_fem_int[vlak],2], method='cubic')
        vz = rbfv_z(exp_coor[vlak,:len(nne),0],exp_coor[vlak,:len(nne),1],exp_coor[vlak,:len(nne),2])

        fem_vlak_rbf[vlak,:len(nne),0] = vx
        fem_vlak_rbf[vlak,:len(nne),1] = vy
        fem_vlak_rbf[vlak,:len(nne),2] = vz

    print(fem_vlak_rbf.shape)
    print("Done!")
    femd.close()
    expd.close()
    # // Below all the information needed is written out:
    # EXP coordinates (x,y and/or z), EXP displacements (x,y and/or z), NUM linear interpolated coordinates (2D/3D),
    #  NUM linear interpolated displacements (2D/3D), NUM RBF interpolated data (displacements (2D/3D), number of EXP nodes,
    # number of NUM nodes, dimension of data (2D,3D)
    # - abbreviations of order above: ce,de,fci,fdi,fvr,nne,nnf,dm
    return(exp_coor,exp_disp,fem_coor_int,fem_disp_int,fem_vlak_rbf,len(nne),len(nn),len(dimension))