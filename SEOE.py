# -----------------------------------------------------------------------------------------------------------------
# This file was written to calculate the Stardard Error of Estimation error for a specific project. The order of the input files are explained
# below. The file can be used for any data if the input files are accoring to the specified format.

# Francè Bresler
# 12 February 2020

# ----------------------------------------------------------------------------------------------------------------

import numpy as np 

# // This function calculates  the Stardard Error of Estimation, takes 2 necessary inputs and 1 optional:
#   -   1. fvr = FEM data after RBF interpolation, contains all the data wanted for calculations.
#           - "fvr" must be 3D: 1 - number of increments, 2 - nodes as rows, 3 - dimensions and other data as columns
#           - "fvr" the first few columns must be the same data as the columns in "de"
#           - in this code is was coded that column 1-2/3 is the x,y,z-disp directions
#   -   2. de = EXP disp data with column 1 as x-direc, 2 as y-direction and 3rd direction optional as z-direction (one direction is also okay)
#   -   3. *args = list of increments to calculate. If args = 0 then the default is to calculate all the increments present in the data
#                       and "fvr" and "de" needs to be 3D arrays. There are 2 input options:
#           - a. only one input is given - int saying the first few increments corresponding to the int will be calculated
#           - b. 2 or more inputs = 1st input is int saying how many increments
#                                   2nd and rest optional more, are the increments to calculate 

def SEOE_disp(fvr,de,*args):
    size = de.shape
    n_dim = size[-1]
    n_nodes = size[-2]
    if len(args) == 0:
        incn = range(size[0])
        levels = size[0]
    elif len(args)==1:
        incn = range(args[0])
        levels = args[0]
    else:
        incn = args[1:]
        levels = args[0]

    SEOE = []
    for d in range(n_dim):
        rsqi = 0
        for j in incn:
            est_sq = 0
            for n in range(n_nodes):
                if np.isnan(fvr[j,n,d]) == True:
                    fvr[j,n,d] = 0.0
                est_sq = est_sq + ((fvr[j,n,d])-(de[j,n,d]))**2
            rqii = np.sqrt(est_sq/(n_nodes - 2))
            rsqi = rsqi + rqii   
        rsq = rsqi/levels
        SEOE.append(rsq)
    # returns all SEOE for the data in one list, the order corresponds to the data in the same order as the columns in "fvr"
    return(SEOE)

# // This function calculates the Stardard Error of Estimation for stress and stretch, taking 2 necessary inputs and 1 optional input:
#   -   1. fvr = FEM data for the stress and stretch only, or the 1st 2 columns need to contain this data at least.
#                3D array, 1- # of increments, 2- nodes as rows, 3- dimensions as columns
#   -   2. ste = EXP data. 3D array, 1- # of increments, 2- nodes as rows, 3- dimensions as columns.
#   -   3. *args = 3 input options:
#       - a. No input =  all the increment data gets processed
#       - b. 1 int input =  only the first few # of increments gets processed
#       - c. 2 or more =  1st int: number of increments to process.
#                         2nd and rest: increments to process, list all of them next.

def SEOE_str(fvr,ste,*args):
    size = ste.shape
    n_dim = size[-1]
    n_nodes = size[-2]
    if len(args) == 0:
        incn = range(size[0])
        levels = size[0]
    elif len(args)==1:
        incn = range(args[0])
        levels = args[0]
    else:
        incn = args[1:]
        levels = args[0]

    SEOE = []
    for d in range(n_dim):
        rsqi = 0
        for j in incn:
            est_sq = 0
            if np.isnan(fvr[j,d]) == True:
                fvr[j,d] = 0.0
            est_sq = est_sq + ((fvr[j,d])-(ste[j,d]))**2
            rqii = np.sqrt(est_sq/(n_nodes - 2))
            rsqi = rsqi + rqii   
        rsq = rsqi/levels
        SEOE.append(rsq)
    # returns all SEOE for the data in one list, the order corresponds to the data in the same order as the columns in "fvr"
    return(SEOE)
