# Wrapper to call DOT from Python

# ------------------------------------------------------------------------------------------------------------------
# This file was edited as needed and explained where necessary.

# Francè Bresler
# 12 February 2020
# ------------------------------------------------------------------------------------------------------------------


import os
import numpy as nm
import ctypes as ct
from ctypes import byref as B
from py_mentat import *
from py_post import *
import subprocess
import pandas as pd
import csv
import matplotlib.pyplot as plt 
import time
from File_paths import filepaths

class dot:
	
	#Set some local constants
	nInfo 	= 0
	nMethod = 0
	nPrint	= 0
	nMinMax	= 0
	nMaxInt	= 20000000
	nmParam = nm.empty(1, float)
	nmRPRM  = nm.zeros(20, float)
	nmIPRM  = nm.zeros(20, int)

	def __init__(self):
		self.dotlib = ct.windll.LoadLibrary("DOT.dll")

	# The DOT wrapper
	def dotcall(self, x, xl, xu, nCons):

		# Reset nInit
		nInit = 0

		#Initailize all array types
		nDvar = x.shape[0]
		ctDVAR	= ct.c_double * nDvar
		ctCONS	= ct.c_double * nCons
		ctRPRM	= ct.c_double * 20
		ctIPRM	= ct.c_int * 20

		#Initialize all arrays
		RPRM = ctRPRM(*(self.nmRPRM))		#Tells dot to use defaults
		IPRM = ctIPRM(*(self.nmIPRM))		#Tells dot to use defaults
		X    = ctDVAR(*(x))				#Initial values
		XL   = ctDVAR(*(xl))			#Lower bounds
		XU   = ctDVAR(*(xu))			#Upper bounds
		G    = ctCONS(*([0.0]*nCons))	#Constraints

		#Initialize constants
		METHOD  = ct.c_int64( self.nMethod )
		NDV     = ct.c_int64( nDvar )
		NCON    = ct.c_int64( nCons )
		IPRINT  = ct.c_int64( self.nPrint )
		MINMAX  = ct.c_int64( self.nMinMax )
		INFO    = ct.c_int64( self.nInfo )
		OBJ     = ct.c_double( 0.0 )
		MAXINT  = ct.c_int64( self.nMaxInt )

		# Call DOT510
		NRWK    = ct.c_int64()
		NRWKMN  = ct.c_int64()
		NRIWD   = ct.c_int64()
		NRWKMX  = ct.c_int64()
		NRIWK   = ct.c_int64()
		NSTORE  = ct.c_int64()
		NGMAX   = ct.c_int64()
		IERR    = ct.c_int64()

		self.dotlib.DOT510(B(NDV), B(NCON), B(METHOD), B(NRWK), B(NRWKMN), B(NRIWD), B(NRWKMX), B(NRIWK), B(NSTORE), B(NGMAX), B(XL), B(XU), B(MAXINT), B(IERR))

		ctRWK	= ct.c_double * NRWKMX.value
		ctIWK	= ct.c_int64 * NRIWK.value
		IWK	= ctIWK( *([0]*NRIWK.value) )
		WK	= ctRWK( *([0.0]*NRWKMX.value) )

		# Call DOT
		# // Here the original dot.py code was edited for personal use. The code can be adjusted as suited to the user.
		# // The iterations and objective lists were created as my own counter and is not necessary. These lists were 
		#	however used in this pipeline.
		itera = 0
		iterations = []
		objective = []
		while (True):
			self.dotlib.DOT(B(INFO),B(METHOD),B(IPRINT), B(NDV),  B(NCON), B(X), B(XL), B(XU), B(OBJ), B(MINMAX), B(G), B(RPRM), B(IPRM), B(WK), B(NRWKMX), B(IWK), B(NRIWK))
			iterations.append(itera)
			objective.append(OBJ.value)
			itera = itera+1
			if ( INFO.value == 0 ) :	# if the optimisation converged, enter loop
				import ast

				# // Open the "iterations" text file to obtain the current design point/ starting point form the list
				#	obtained by the LHC function
				filec = open("iterations.txt","r")  
				it = int(filec.readline())
				filec.close
				# //
				path1 = filepaths("path1")

				# // Read what the original starting point was of this optimisation run.
				filen = "starting_points.txt"
				filp = os.path.join(path1,filen)
				start = open(filp,'r')
				xxx = ast.literal_eval(start.readlines()[it])
				start.close()

				xxx = nm.array(xxx)
				xc = nm.array(X)
				print(xxx)
				print(xc)
				xa = xc*xxx 	# the final design point from dot is multiplied with the starting point to obtain the 
				print(xa)		# material coefficient values, since the current values are the unbiased values.

				from functions import append_val, expdata, fem_orig_data, final_points	# call the output functions needed
				if len(X)==3:
					xf = [xa[0], xa[1], xa[2]]	# store the optimised point in a form easily written to a text file
				elif len(X)==2:
					xf = [xa[0], xa[1]]
				time.sleep(1)
				append_val(xf,iterations,objective)	# Writes out the iteration file to show how it converged

				final_points(xf)	# store the optimised point in a separate file
				time.sleep(1)

				from functions import material2d, material3d#, material3d_ogden
				if len(X)==3:
					material3d(xa)	# Create new procedure file for the optimised point
				elif len(X)==2:
					material2d(xa)

				time.sleep(1)
				filem = "mat.proc"
				from code_settings import pipeline_files
				pf = pipeline_files()
				p = subprocess.Popen([filepaths(pf.mentat_version),filem], bufsize=2048)  # Start MSC Marc and load the procedure 
				# file which will open the correct NUMERICAL model and change the material properties, start Marc 
				# solver and to save the post file for the current DOT increment, close Marc and continue with the 
				# code below
				p.wait()
				time.sleep(5)

				# Ensure that Marc file converged
				sts = filepaths("fem_sts")
				conver = pd.read_csv(sts, header=None, sep=' ', names=list(range(11)), keep_default_na=False)    
				# open sts file of the optimisation NUMERICAL model, it contains the exit code from Marc.
				
				time.sleep(1)
				c = int(conver.iloc[-3,-1])
				if c == 3004:
					g = -1       # 3004 says the FEM converged and the constraint is satisfied
				else:
					g = 1        # Any other exit number says the FEM did not converge and therefore the constraint 
					# was not satisfied. This is a fail save to ensure the optimised point does adhere to the constraints.
					xv = xa

					time.sleep(1)
					from functions import violated_constr
					violated_constr(xv,c)     # This function saves the parameters which caused non-convergence and 
					# also what the exit number was
				print(g)

				from RBF import RBF_int		# call the RBF function
				fname1,fname2,fname3,fname4,fnamef,fnamee = filepaths("fem_out","pointfem","exp_out","pointexp","fem_dat","exp_dat")

				# // All the data after interpolation
				ce,de,fci,fdi,fvr,nne,nnf,dm = RBF_int(fname1,fname2,fname3,fname4,fnamef,fnamee,g)
				time.sleep(1)
				# //

				expdata(ce,de,fvr,nne)		# Store the experimental data for the starting point
				fem_orig_data(fci,fdi,nnf,nne)	# Store the ouput data for the optimised point's simulation
				time.sleep(1)
				from functions import objectfunc
				objectfunc(objective,iterations)
				
				rem = filepaths("mat_proc_path")
				os.remove(rem)
				# os.remove('../Cylinder/mat.proc')
				break
			else: # if the optimisation procedure haven't converged yet
				fname1 = filepaths("fem_out")    
				os.remove(fname1)     #delete the current DOT optimisation increment's NUMERICAL data
				self.evaluate(X, OBJ, G, self.nmParam)
		# //

		rslt = nm.empty( 2+nDvar, float)
		rslt[0] = OBJ.value
		rslt[1] = 0.0
		if len(G) > 0 :
			rslt[1] = max(G)
		for i in range( nDvar ):
			rslt[2+i] = X[i]
		return rslt

	def evaluate(self, x, obj, g, param):
		obj.value = 2.0*(x[0]*x[1] + x[0]*x[2] + 2.0*x[1]*x[2])
		g[0] = 1.0 - 0.5*x[0]*x[1]*x[2]
		return
