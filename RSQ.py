# -----------------------------------------------------------------------------------------------------------------
# This file was written to calculate the R^2 error for a specific project. The order of the input files are explained
# below. The file can be used for any data if the input files are accoring to the specified format.

# Francè Bresler
# 12 February 2020

# ----------------------------------------------------------------------------------------------------------------

import numpy as np 

# // This function calculates the R^2 of the data for the dispalacements, taking 2 necessary inputs and 1 optional input:
#   -   1. fvr = FEM data for the displacements. 3D array, 1- # of increments, 2- nodes as rows, 3- dimensions as columns
#   -   2. de = EXP data. 3D array, 1- # of increments, 2- nodes as rows, 3- dimensions as columns.
#   -   3. *args = 3 input options:
#       - a. No input =  all the increment data gets processed
#       - b. 1 int input =  only the first few # of increments gets processed
#       - c. 2 or more =  1st int: number of increments to process.
#                         2nd and rest: increments to process, list all of them next.
def RSQ_disp(fvr,de,*args):
    size = fvr.shape
    dm = de.shape[-1]
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

    RSQ = []
    alp = []
    bet = []
    for d in range(dm):
        rsq = 0
        sum_xy = 0
        sum_x = 0
        sum_y = 0
        sum_x2 = 0
        alpha = 0
        beta = 0
        for j in incn:
            dx = de[j,:,d]
            dy = fvr[j,:,d]
            sum_x = np.sum(dx)
            sum_y = np.sum(dy)
            sum_xy = np.sum((dx*dy))
            sum_x2 = np.sum((dx**2))
            alpha = ((n_nodes*sum_xy)-(sum_x*sum_y))/((n_nodes*sum_x2)-(sum_x)**2)
            beta = (sum_y-(alpha*sum_x))/n_nodes
            dy_est = (alpha*dx) + beta
            mean_y = np.mean(dy)
            est_sq = (dy_est-mean_y)**2
            act_sq = (dy-mean_y)**2
            rsqi = sum(est_sq)/sum(act_sq)
            if np.isnan(rsqi) == True:
                rsqi = 0
            rsq = rsq + rsqi
        rsq = rsq/levels
        alp.append(alpha)
        bet.append(beta)
        RSQ.append(rsq)
    # returns 3 lists   -  line of best fit eq : y = alpha*x + beta 
    # 1 = list of the R^2 of each dimension
    # 2 = list of the alpha values in the linear equation
    # 3 = list of the beta values in the linear equation
    return(RSQ,alp,bet)

# // This function calculates the Root Mean Square of the data for the stress and strain, taking 2 necessary inputs and 1 optional input:
#   -   1. fvr = FEM data for the stress and stretch only, or the 1st 2 columns need to contain this data at least.
#                3D array, 1- # of increments, 2- nodes as rows, 3- dimensions as columns
#   -   2. ste = EXP data. 3D array, 1- # of increments, 2- nodes as rows, 3- dimensions as columns.
#   -   3. *args = 3 input options:
#       - a. No input =  all the increment data gets processed
#       - b. 1 int input =  only the first few # of increments gets processed
#       - c. 2 or more =  1st int: number of increments to process.
#                         2nd and rest: increments to process, list all of them next.
def RSQ_st(fvr,ste,*args):
    size = fvr.shape
    n_dim = ste.shape[-1]
    n_nodes = size[0]
    if len(args) == 0:
        levels = size[-2]
    elif len(args)==1:
        levels = args[0]
    else:
        inc = args[1:]
        levels = args[0]
        dxi = np.zeros([levels,n_dim])
        dyi = np.zeros([levels,n_dim])
        for j in range(len(inc)):
            dxi[j,:] = ste[inc[j],:]
            dyi[j,:] = fvr[inc[j],:]
        
    RSQ = []
    alp = []
    bet = []
    for d in range(n_dim):
        rsq = 0
        sum_xy = 0
        sum_x = 0
        sum_y = 0
        sum_x2 = 0
        alpha = 0
        beta = 0

        if len(args)<2:
            dx = ste[:levels,d]
            dy = fvr[:levels,d]
        else:
            dx = dxi[:,d]
            dy = dyi[:,d]
        sum_x = np.sum(dx)
        sum_y = np.sum(dy)
        xy = dx*dy
        x2 = dx**2
        sum_xy = np.sum(xy)
        sum_x2 = np.sum(x2)
        alpha = ((n_nodes*sum_xy)-(sum_x*sum_y))/((n_nodes*sum_x2)-(sum_x)**2)
        beta = (sum_y-(alpha*sum_x))/n_nodes
        dy_est = (alpha*dx) + beta
        mean_y = np.mean(dy)
        est_sq = (dy_est-mean_y)**2
        act_sq = (dy-mean_y)**2

        rsqi = sum(est_sq)/sum(act_sq)
        if np.isnan(rsqi) == True:
            rsqi = 0
        alp.append(alpha)
        bet.append(beta)
        RSQ.append(rsqi)
    # returns 3 lists   -  line of best fit eq : y = alpha*x + beta 
    # 1 = list of the R^2 of each dimension
    # 2 = list of the alpha values in the linear equation
    # 3 = list of the beta values in the linear equation
    return(RSQ,alp,bet)