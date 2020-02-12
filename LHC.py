# -----------------------------------------------------------------------------------------------------------------
# This file uses the pyDOE module to use the lhc function and obtain 10 random starting points for the optimisation
# procedure. Since random points between 0-1 are chosen, these points need to be scaled between the upper and lower 
# bound given.

# These upper and lower bounds can be changed as desired, by default they were chosen to be 20% above and below the 
# original material model's values used in the "Experimental" model:
#       orig_coef = [0.26056762, 0.09754981, 0.05750069]

# Francè Bresler
# 12 February 2020
# ---------------------------------------------------------------------------------------------------------------
from pyDOE import *
import numpy as np 

# // The latent hyper cube (LHC) function creates "s" number of starting points of dimension "n", these starting points are chosen
# equally in the design space between 0 and 1.
#   - Since the coefficients are not values between 0 and 1, lower and upper bounds are specified and the values onbtained
#     from the LHC function, is scaled between the specified lower and upper bounds.
def LHC(n,s):
    cl = [0.208454096,0.078039848,0.046000552] # lower bounds for the material coefficients
    cu = [0.312681144,0.117059772,0.069000828] # upper bounds for the material coefficients
    
    samp = lhs(n, samples=s, criterion="maximin")

    # the for loop scales the values to be equally spaced between the upper and lower bounds specified above
    for i in range(s):
        for j in range(n):
            samp[i,j] = samp[i,j]*(cu[j]-cl[j]) + cl[j]

    return samp,cl,cu
