# -----------------------------------------------------------------------------------------------------------------
# Using the RMS value to determine the correlation between the NUM data and the EXP data,
# The final RMS value is used as the objective function during the optimisation.

# Francè Bresler
# 12 February 2020

# -----------------------------------------------------------------------------------------------------------------



import numpy as np 

# //RMS_disp takes 2 necessary inputs and 1 optional input:
#   1. exp_data = the EXP data for the displacements output from the RBF function (must only have x,y and/or z columns)
#   2. fem_data = the NUM data from the RBF function with all the interpolated values (namely the 4/5 column output after RBF interp.)
#   3. *args = this input is optional with the following choices:
#       -   no input: the exp_data and fem_data arrays need to be 3D, the function will work out the RMS for all the simulation increments
#           ex. RMS_disp([5,10,2/3],[5,10,2/3 (only disp) or 4/5(disp and stress/strain)])
#               all 5 increments will be used to calculate 1 RMS value
#       -   1 input: the first few increments corresponding to the intreger value
#       -   2 or more values: Value 1 = intreger value saying how many increments (1 or more)
#                             Value 2 and more = the increments to use
#           ex. RMS(exp_data,fem_data,1,15) = only increment 15 will be calculated
#           ex. RMS(exp_data,fem_data,3,8,10,12) = increment 8, 10 and 12 are the only 3 increments used
def RMS_disp(exp_data,fem_data,*args):
    size = exp_data.shape
    n_inc = size[0]
    n_nodes = size[-2]
    n_dim = size[-1]
    if len(args) == 0:
        n_inc = range(size[0])
        levels = size[0]
    elif len(args) == 1:
        n_inc = range(args[0])
        levels = args[0]
    else:
        n_inc = args[1:]
        levels = args[0]

    rmsi = 0
    for i in n_inc:
        rmsn = 0
        for k in range(n_dim):
            node_inc = n_nodes
            diff = 0
            for j in range(node_inc):
                df = fem_data[i,j,k]
                de = exp_data[i,j,k]
                if np.isnan(df) == True:        #Make sure non of the data contains NaN values
                    node_inc = node_inc-1
                    continue
                else:
                    diff = diff + (df-de)**2    # Sum all the disp differences together for the dimension
            # avg = abs(np.average(exp_data[i,:,k]))
            avg = max(abs(exp_data[i,:,k]))
            if avg == 0.0:
                avg = 1.0
            rmsid = np.sqrt((diff/node_inc))    # calc the RMS for the direction in the specific increment
            nom = rmsid/avg
            rmsn = rmsn + nom          # RMS for the increment is equal to summation of each normalised direction RMS
        rmsi = rmsi + rmsn        # Sum all the increment RMS's 
    rms = rmsi/levels          # Final Objection function RMS = avg of all the increment RMS's

    # calculate only each direction's RMS
    rmsd = []
    for k in range(n_dim):
        rmsi = 0
        rmsn = 0
        for i in n_inc:
            node_inc = n_nodes
            diff = 0
            for j in range(node_inc):
                df = fem_data[i,j,k]
                de = exp_data[i,j,k]
                if np.isnan(df) == True:        # make sure there is no NaN values
                    node_inc = node_inc-1
                    continue
                else:
                    diff = diff + (df-de)**2    # Sum all the disp differences together
            # avg = abs(np.average(exp_data[i,:,k]))
            avg = max(abs(exp_data[i,:,k]))
            if avg == 0.0:
                avg = 1.0
            rmsid = np.sqrt((diff/node_inc))    # Calculate the increment RMS for the direction
            nom = rmsid/avg
            rmsn = rmsn + nom               # Sum each increment's RMS  
        rmsi = rmsn/levels      # Direction's RMS value gets averaged over all the increments
        rmsd.append(rmsi)       # Save the direction's RMS
    
    return(rms,rmsd) # - rms = total displacement RMS
                     # - rmsd = each displacement direction's RMS value

# //RMS_stvssr takes 2 necessary inputs and 1 optional input:
#   1. exp_data = the EXP data for the stress and/or stretch output from the RBF function 
#                   (must only have stress and/or stretch column NOT the nodes column as well)
#   2. fem_data = the NUM data from the RBF function with all the interpolated values 
#                   (must only have stress and/or stretch column NOT the displacements columns as well)
#   3. *args = this input is optional with the following choices:
#       -   no input: the exp_data and fem_data arrays need to be 3D, the function will work out the RMS for all the simulation increments
#           ex. RMS_stvssr([5,10,1/2],[5,10,1/2])
#               all 5 increments will be used to calculate 1 RMS value for stress or 2 RMS values one for stress & 1 for strain
#       -   2 or more values: Value 1 = int value saying how many increments (1 or more)
#                             Value 2 and more = the increments to use
#           ex. RMS(exp_data,fem_data,1,15) = only increment 15 will be calculated
#           ex. RMS(exp_data,fem_data,3,8,10,12) = increment 8, 10 and 12 are the only 3 increments used
def RMS_stvssr(exp_data,fem_data,*args):
    size = fem_data.shape
    n_inc = size[0]
    dm = size[-1]
    if len(args) == 0:
        n_nodes = range(size[-2])
        levels = size[-2]
    elif len(args) == 1:
        n_nodes = range(args[0])
        levels = args[0]
    else:
        n_nodes = args[1:]
        levels = args[0]

    rms = []
    for d in range(dm):
        diff = 0
        for i in n_nodes:
            node_inc = n_inc
            df = fem_data[i,d]
            de = exp_data[i,d]
            if np.isnan(df) == True:        # Make sure there are no NaN values
                node_inc = node_inc-1
                continue
            else:
                diff = diff + (df-de)**2      # Sum all the Vector value differences together
        # avg = abs(np.average(exp_data[:,d]))
        avg = max(abs(exp_data[:,d]))
        if avg == 0.0:
            avg = 1.0
        rmsie = np.sqrt((diff/node_inc))            # Calculate the Vector RMS for the increment
        nom = rmsie/avg                       # Normalise the increment RMS with the avg
        rmsst = (nom/levels)               # Obtain the final RMS value as the avg of all the increments
        rms.append(rmsst)

    return(rms) 



    