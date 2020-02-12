# -----------------------------------------------------------------------------------------------------------------
# This file serve as an example file and is not completely necessary.
# Three graphs are drawn here after the optimisation procedure, with a 4th one as optional if two algorihms are tested
#
# In the following lines it is explaned how and where the code needs to change for 2 algorithms:
#       - lines 52 - 53
#       - lines 228 - 237
#       - lines 486 - 500
#
# Francè Bresler
# 12 February 2020
#
# ------------------------------------------------------------------------------------------------------------------

# All necessary functions and files

from py_post import *
from py_mentat import *
import numpy as np 
import pandas as pd 
import os
import matplotlib.pyplot as plt
import matplotlib.lines as pll
import csv
import time
import ast
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

from functions import append_val, final_points
from SEOE import SEOE_disp, SEOE_str
from RMS import RMS_disp, RMS_stvssr
from RSQ import RSQ_disp, RSQ_st

from File_paths import filepaths
 
# // This file serve as an example how the rest of the pipeline and saved data was used for usefull graphs. 
#   - The engineering stress was determined by using the Mooney-Rivlin three parameter model for a uni-axial tensile case.
#   - A stretch range of 0.4 - 3.0 stretch was used to determine the engineering stress for the MR model. This range
#     was chosen to determine the predictability of the material model parameters obtained by the optimisation 
#     for compression and tension, for the specified silicone rubber used.
#   - The "EXP" model used the material parameters in the "orig_mat_param" variable, which was obtained from a previous
#     study on silicone rubber.
def post_process(dm,*args):
    print('Start')

#---------------------------------------------------------------------------------------
# ---- Unindent the following lines when two optimisation algorithms are used

# NB !!!!! - Remember to unindent lines 

    # algpat1 = "SQP Optimisation Algorithm"
    # algpat2 = "SLP Optimisation Algorithm"
#---------------------------------------------------------------------------------------
    
    # plt.close('all')
    filen = "final_points.txt" 
    if len(args)==0:
        path1 = filepaths("path1")
        path2 = filepaths("path2")

        filet = open("iterations.txt","r")  # determine how many starting point iterations there were (ex. 10), 
                                            # The results for each starting point iteration will be plotted. 
        it = int(filet.readline())
        filet.close
    elif len(args)>0:
        pat = filepaths("path1")
        pat2 = filepaths("path2")
        fl = args[0]
        path1 = os.path.join(pat,fl)
        path2 = os.path.join(pat2,fl)

        it = 0

    filepath = os.path.join(path1,filen)
    if not os.path.exists(path1):
        os.makedirs(path1)

    filepath2 = os.path.join(path2,filen)
    if not os.path.exists(path2):
        os.makedirs(path2)

    filere = "Final_results.txt"            # Create a text file with the final errors 
    filepathf = os.path.join(path1,filere)  
    filerb = "Best_Final_results.txt"       # Create a text file with only the optimisation run with
                                            # the best smallest objective function
    filepathb = os.path.join(path1,filerb)

    file = open(filepath,'r')
    s = ast.literal_eval(file.readlines()[0])
    xf = np.empty([it+1,len(s)])
    rms_d = []

    for i in range(it+1):
        file.seek(0,os.SEEK_SET)
        x = ast.literal_eval(file.readlines()[i])   # read from the "final_points.txt" read the optimised parameters
        for j in range(len(x)):
            xf[i,j] = x[j]

    coefshape = xf.shape[1]
    
    orig_coef = [0.26056762, 0.09754981, 0.05750069]    # The original parameters as used by the "EXP" model
    lower_coef = [0.208454096,0.078039848,0.046000552]  # Lower bound
    high_coef = [0.312681144,0.117059772,0.069000828]   # Upper bound

    ocdl = []
    ocdh = []

    # // Determine the difference in the coefficients and bounds for the % calculation
    for hl in range(len(orig_coef)):
        ocdl.append(orig_coef[hl]-lower_coef[hl])
        ocdh.append(high_coef[hl]-orig_coef[hl])
    per_acc = np.empty([it+1,len(orig_coef)])
    # //

    # // Using the Mooney-Rivlin three parameter equation for a uni-axial tensile case, the engineering stress is 
    # calculated for a given stretch range
    expsr_e = np.linspace(0.4,3.0,106)          # Given stretch range
    expsr_f = np.linspace(0.4,3.0,106)
    stef = np.empty([it+1,len(expsr_e),2])
    fvrsf = np.empty([it+1,len(expsr_e),2])
    expse_f = np.empty([it+1,len(expsr_e),2])    # Empty arrays to store the calculated eng stresses
    expsf_f = np.empty([it+1,len(expsr_e),2])
    expsff_f = np.empty([it+1,len(expsr_e),2])

    rms_stsr = []
    for k in range(it+1):               # Loop to calculate each optimisation run's results, errors and graphs
        fname1 = "RBF_data%s.npy"%k     # Import the saved "NUM" results, after interpolation
        fname2 = "exp_data%s.npy"%k     # Import the "EXP" results
        filepath1 = os.path.join(path1,fname1)
        filepath2 = os.path.join(path1,fname2)
        RBFdata = np.load(filepath1)
        expdata = np.load(filepath2)
        der = expdata[:,:,dm:(dm*2)]    # The first few columns are the coordinates, the next few are the displacements
                                        # which are accessed here
        obj = RMS_disp(der,RBFdata)     # The RMS error in the displacements are calculated
        rms_d.append(obj[0])

        # // Using the MR three parameter equation for uni-axial tensile case, 
        # - the eng stress for the "EXP" model is calculated = expst_e
        # - the eng stress for the "NUM" model is calculated = expst_f
        ty = 'Mooney-Rivlin'
        expst_e = 2*orig_coef[0]*(expsr_e-(1/(expsr_e**2))) + 4*orig_coef[2]*((expsr_e**2)+(2/expsr_e)-3)*(expsr_e-(1/(expsr_e**2))) + 2*orig_coef[1]*(1-(1/(expsr_e**3)))
        expst_f = 2*xf[k,0]*(expsr_f-(1/(expsr_f**2))) + 4*xf[k,2]*((expsr_f**2)+(2/expsr_f)-3)*(expsr_f-(1/(expsr_f**2))) + 2*xf[k,1]*(1-(1/(expsr_f**3)))

        # The stretch values in the "NUM" model is calculated using the eng stress results above for the "EXP" model
        expsr_ff = np.zeros(len(expst_e))
        for i in range(len(expst_e)):
            xc = lambda x: expst_e[i] - (2*xf[k,0]*(x-(1/(x**2))) + 4*xf[k,2]*((x**2)+(2/x)-3)*(x-(1/(x**2))) + 2*xf[k,1]*(1-(1/(x**3))))
            sol = fsolve(xc,expsr_e[i])
            # print(sol)
            expsr_ff[i] = sol
        # //

        ste = np.zeros([len(expsr_e),2])
        fvrs = np.zeros([len(expsr_e),2])

        # // Storing each optimisation run's eng stress and stretch values
        expse_f[k,:,0], expse_f[k,:,1] = expst_e, expsr_e   # Storing the "EXP" model's eng stress with the given stretch range
        expsf_f[k,:,0], expsf_f[k,:,1] = expst_f, expsr_f   # Storing the "NUM" model's eng stress for the same given stretch range
        expsff_f[k,:,0], expsff_f[k,:,1] = expst_f, expsr_ff    # Storing the "NUM" model's eng stress, with the stretch values calculated
                                                                # at the "EXP" eng stress values, to calculate the error between the stretch results
        
        # // Storing only the current optimisation run's eng stress and stretch values for the error calculations
        ste[:,0],ste[:,1] = expst_e, expsr_e
        fvrs[:,0],fvrs[:,1] = expst_f, expsr_ff
        stef[k,:,:] = ste
        fvrsf[k,:,:] = fvrs
        rmsr = RMS_stvssr(ste,fvrs) # determine the RMS error for the stress and then the stretch
        rms_stsr.append(sum(rmsr))  # sum the two RMS errors for a combined RMS error for the stress and stretch

        # // Determine the % how far the optimised parameter is from the original parameter, 
        # with the bounds being 0 % and the orig_coef parameter 100 %
        for l in range(len(orig_coef)):
            if (xf[k,l] <= orig_coef[l]):
                diff = orig_coef[l] - xf[k,l]
                per = (abs(diff - ocdl[l])/ocdl[l])*100
                per_acc[k,l] = per
            elif (xf[k,l] > orig_coef[l]):
                diff = xf[k,l] - orig_coef[l]
                per = (abs(diff - ocdh[l])/ocdh[l])*100
                per_acc[k,l] = per 
        # //          

    # // Determine which optimisation run had the best RMS value for the displacements and set that run as the preliminary "best",
    # for the optimisation runs which have RMS values within 5 % of the best, the final best optimisation run will be 
    # chosen for the run which have the best combined stress and stretch RMS value.  
    b_rms = min(rms_d)
    bestd = np.where(np.array(rms_d) == b_rms)[0][0]
    ub = rms_d[bestd]*0.05 + rms_d[bestd]
    lb = rms_d[bestd]*0.05 - rms_d[bestd]
    cnt = []
    rms_obj = []
    for t in range(it+1):
        if (rms_d[t]>=lb) and (rms_d[t]<=ub):
            cnt.append(t)
    if len(cnt)>1:
        for tt in cnt:
            rms_obj.append(rms_stsr[tt])
        rb_rms = min(rms_obj)
        be = np.where(np.array(rms_obj) == rb_rms)[0][0]
        best = cnt[be]
    else:
        best = best
    # //

    print(best)
    fn = 2
    objl = []
    objlSP = []

    # // Generate the graphs and error results for each optimisation run
    for bb in range(it+1):
        fname = "fem_orig_data%s.npy"%bb
        fname1 = "RBF_data%s.npy"%bb
        fname2 = "exp_data%s.npy"%bb            # Import all the saved data

        filepath = os.path.join(path1,fname)
        filepath1 = os.path.join(path1,fname1)
        filepath2 = os.path.join(path1,fname2)  # obtain the correct file paths

# -----------------------------------------------------------------------------------------------------------------
    #----- Unindent the following lines only if two optimisation algorithms are used
    # ---- Only unindent before the 2nd algorithm is started, NOT during the 1st algorithm run
    
    # NB !!!! --- Remember to unindent lines 486 - 500
        
        # fileobj = "objectfunction%s.npy"%bb
        # fileiter = "niter%s.npy"%bb
        # patobj = os.path.join(path1,fileobj)
        # patiter = os.path.join(path1,fileiter)
        # patobjSP = os.path.join(path2,fileobj)
        # objectivefunc = np.load(patobj)
        # objl.append(objectivefunc[-1])
        # objectivefuncSP = np.load(patobjSP)
        # objlSP.append(objectivefuncSP[-1])
        # numiter = np.load(patiter)
# -----------------------------------------------------------------------------------------------------------------
        femdata = np.load(filepath)
        RBFdata = np.load(filepath1)
        expdata = np.load(filepath2)

        sf = RBFdata.shape[0]
        fvr = RBFdata
        de = expdata[:,:,dm:(dm*2)]     # obtain only the "EXP" model's disp data and not coordinate data
        seoefd = SEOE_disp(fvr,de)

        ste = stef[bb,:,:]
        fvrs = fvrsf[bb,:,:]
        expsr_e = expse_f[bb,:,1]
        expst_e = expse_f[bb,:,0]      # obtain the engineering stress and stretch calculated
        expsr_f = expsf_f[bb,:,1]
        expst_f = expsf_f[bb,:,0]

        rmsr = RMS_stvssr(ste,fvrs)         # determine the RMS error for eng stress and stretch
        rsqs,alps,bets = RSQ_st(fvrs,ste)   # determine the R^2 error for eng stress and stretch
        # print(rsqs)
        # print(alps)
        # print(bets)
        seoef = SEOE_str(fvrs,ste)      # determine the Standard Error of Estimation for eng stress and stretch
        incn = sf
        rmsd = RMS_disp(de,fvr)     # determine the RMS error for each displacement direction
        rsq = RSQ_disp(fvr,de)      # determine the R^2 error for each displacement direction

        # // Error result summary are written to text files
        if bb == best:
            with open(filepathb,"a+",encoding="utf-8") as fr:
                fr.write("According to THE BEST POINT - "+ str(bb) +" the results are: \n")
                fr.write("The optimal %s coefficients are: \n"%ty)
                fr.write("\n")
                if len(orig_coef) == 3:
                    fr.write("C10 = " + str(xf[bb,0]) + "\nC01 = " + str(xf[bb,1]) + "\nC20 = " + str(xf[bb,2]))
                else:
                    # web for unicodes for greek letter = https://pythonforundergradengineers.com/unicode-characters-in-python.html
                    fr.write('\u03BC1 = ' + str(xf[bb,0]) + "\n\u03BC2 = " + str(xf[bb,1]) + "\n\u03BC3 = " + str(xf[bb,2]) + '\n\u03B11 = ' + str(xf[bb,3]) + "\n\u03B12 = " + str(xf[bb,4]) + "\n\u03B13 = " + str(xf[bb,5]))
                fr.write("\n")
                fr.write("Percentage accuracy towards original model's coefficients: \n")
                fr.write("\n")
                if len(orig_coef)==6:
                    fr.write("\u03BC1 = " + str(per_acc[bb,0]) + "\n\u03BC2 = " + str(per_acc[bb,1]) + "\n\u03BC3 = " + str(per_acc[bb,2]) + "\n\u03B11 = " + str(per_acc[bb,3]) + "\n\u03B12 = " + str(per_acc[bb,4]) + "\n\u03B13 = " + str(per_acc[bb,5]))
                if len(orig_coef)==3:
                    fr.write("C10 = " + str(per_acc[bb,0]) + "\nC01 = " + str(per_acc[bb,1]) + "\nC20 = " + str(per_acc[bb,2]))
                fr.write("\n")
                fr.write("It produced  the following errors: \n")
                fr.write("\n")
                fr.write("Displacement RMS = " + str(rmsd[0]) + "\n")
                fr.write("Stress RMS = " + str(rmsr[0])+ "\n")
                fr.write("Strain RMS = " + str(rmsr[1])+ "\n")
                fr.write("Displacement x R**2 = " + str(rsq[0][0]) + "\n")
                fr.write("Displacement y R**2 = " + str(rsq[0][1]) + "\n")
                fr.write("Displacement z R**2 = " + str(rsq[0][2]) + "\n")
                fr.write("Stress R**2 = " + str(rsqs[0]) + "\n")
                fr.write("Strain R**2 = " + str(rsqs[1]) + "\n")
                fr.write("Displacement x Standard Error of Estimate = " + str(seoefd[0]) + "\n")
                fr.write("Displacement y Standard Error of Estimate = " + str(seoefd[1]) + "\n")
                fr.write("Displacement z Standard Error of Estimate = " + str(seoefd[2]) + "\n")
                fr.write("Stress Standard Error of Estimate = " + str(seoef[0]) + "\n")
                fr.write("Strain Standard Error of Estimate = " + str(seoef[1]) + "\n")
                fr.write("\n")
                fr.close()
        else:
            with open(filepathf,"a+",encoding="utf-8") as fr:
                fr.write("According to Starting point - "+ str(bb) +" the results are: \n")
                fr.write("The optimal %s coefficients are: \n"%ty)
                fr.write("\n")
                if len(orig_coef) == 3:
                    fr.write("C10 = " + str(xf[bb,0]) + "\nC01 = " + str(xf[bb,1]) + "\nC20 = " + str(xf[bb,2]))
                else:
                    # web for unicodes for greek letter = https://pythonforundergradengineers.com/unicode-characters-in-python.html
                    fr.write('\u03BC1 = ' + str(xf[bb,0]) + "\n\u03BC2 = " + str(xf[bb,1]) + "\n\u03BC3 = " + str(xf[bb,2]) + '\n\u03B11 = ' + str(xf[bb,3]) + "\n\u03B12 = " + str(xf[bb,4]) + "\n\u03B13 = " + str(xf[bb,5]))
                fr.write("\n")
                fr.write("\nPercentage accuracy towards original model's coefficients: \n")
                fr.write("\n")
                if len(orig_coef)==6:
                    fr.write("\u03BC1 = " + str(per_acc[bb,0]) + "\n\u03BC2 = " + str(per_acc[bb,1]) + "\n\u03BC3 = " + str(per_acc[bb,2]) + "\n\u03B11 = " + str(per_acc[bb,3]) + "\n\u03B12 = " + str(per_acc[bb,4]) + "\n\u03B13 = " + str(per_acc[bb,5]))
                if len(orig_coef)==3:
                    fr.write("C10 = " + str(per_acc[bb,0]) + "\nC01 = " + str(per_acc[bb,1]) + "\nC20 = " + str(per_acc[bb,2]))
                fr.write("\n")
                fr.write("\nIt produced  the following errors: \n")
                fr.write("\n")
                fr.write("Displacement RMS = " + str(rmsd[0]) + "\n")
                fr.write("Stress RMS = " + str(rmsr[0])+ "\n")
                fr.write("Strain RMS = " + str(rmsr[1])+ "\n")
                fr.write("Displacement x R**2 = " + str(rsq[0][0]) + "\n")
                fr.write("Displacement y R**2 = " + str(rsq[0][1]) + "\n")
                fr.write("Displacement z R**2 = " + str(rsq[0][2]) + "\n")
                fr.write("Stress R**2 = " + str(rsqs[0]) + "\n")
                fr.write("Strain R**2 = " + str(rsqs[1]) + "\n")
                fr.write("Displacement x Standard Error of Estimate = " + str(seoefd[0]) + "\n")
                fr.write("Displacement y Standard Error of Estimate = " + str(seoefd[1]) + "\n")
                fr.write("Displacement z Standard Error of Estimate = " + str(seoefd[2]) + "\n")
                fr.write("Stress Standard Error of Estimate = " + str(seoef[0]) + "\n")
                fr.write("Strain Standard Error of Estimate = " + str(seoef[1]) + "\n")
                fr.write("\n")
                fr.close()
            # //

        # // the final simulation increment is used for the displacement data, this is to obtain the final results 
        # at full deformation  
        inc = incn-1

        x_new = expdata[inc,:,0]
        y_new = expdata[inc,:,1]
        dx_new = expdata[inc,:,3]
        dy_new = expdata[inc,:,4]
        z_new = expdata[inc,:,2]
        dz_new = expdata[inc,:,5]
        x_rbf = RBFdata[inc,:,0]
        y_rbf = RBFdata[inc,:,1]
        z_rbf = RBFdata[inc,:,2]
        orgx = expdata[0,:,0]
        orgy = expdata[0,:,1]
        px = orgx + x_rbf
        py = orgy + y_rbf
        # //

        # // The eng stress and stretch is determined for the increment
        x_strain_exp = np.linspace(0.4,3.0,50)
        y_stress_exp = 2*orig_coef[0]*(x_strain_exp-(1/(x_strain_exp**2))) + 4*orig_coef[2]*((x_strain_exp**2)+(2/x_strain_exp)-3)*(x_strain_exp-(1/(x_strain_exp**2))) + 2*orig_coef[1]*(1-(1/(x_strain_exp**3)))
        y_stress_fem = 2*xf[bb,0]*(x_strain_exp-(1/(x_strain_exp**2))) + 4*xf[bb,2]*((x_strain_exp**2)+(2/x_strain_exp)-3)*(x_strain_exp-(1/(x_strain_exp**2))) + 2*xf[bb,1]*(1-(1/(x_strain_exp**3)))
        # //

        # // The error for the eng stress and stretch in the increment
        rsqid,alpd,betd = RSQ_disp(fvr,de,1,inc)
        rms_f,rmsd_f = RMS_disp(de,fvr,1,inc) 
        seoedd = SEOE_disp(fvr,de,1,inc)
        # //

        if bb == best:
            bb1 = bb
            bb = "best"
        else:
            bb = bb
            bb1 = bb
        # Line2D([0], [0], marker='o', color='b', label='FEM-data', markerfacecolor='w', markersize=3)

        # // The total error for the eng stress and stretch
        rmstot = RMS_stvssr(ste[:,:],fvrs[:,:])
        rsqtot,alptot,bettot = RSQ_st(fvrs[:,:],ste[:,:])
        seoetot = SEOE_str(fvrs[:,:],ste[:,:])
        # //

        # // Figure 1 plots the Eng stress vs. stretch for the optimisation run
        plt.figure(fn, figsize=(11.0,8.5))
        plt.plot(expsr_e,expst_e,"r-")
        plt.plot(expsr_f,expst_f,"b-")
        yerr = np.std(y_stress_exp)/np.sqrt(len(y_stress_exp))
        plt.plot(x_strain_exp,y_stress_exp,'g-')
        plt.errorbar(x_strain_exp,y_stress_exp, yerr=yerr, fmt='none', ecolor='g', elinewidth=0.5, capsize=2, errorevery=2)
        plt.xlabel('Stretch',fontsize=18)
        plt.ylabel('Engineering Stress [MPa]',fontsize=18)
        plt.legend(['EXP Model Parameters','NUM Model Parameters','EXP Model Actual Stretch Range'],fontsize=14)
        plt.text(0.7*max(expsr_e), max(expst_e)-0.7*(abs(max(expst_e)-min(expst_e))), "$R^2$ = %s, SEOE = %s\n RMS = %s"%(round(rsqtot[0],4),round(seoetot[0],4),round(rmstot[0],4)))
        plt.grid("on")
        fig5 = "Figure1_%s.pdf"%bb
        figpath5 = os.path.join(path1,fig5)
        plt.savefig(figpath5, dpi=600)
        plt.close(fig=fn)
        fn += 1

        # // Figure 2 plots the x, y and z displacement data to the coordinates
        plt.figure(fn, figsize=(15.0,8.5))
        plt.figure(fn).legend(handles=[pll.Line2D([0], [0], marker='o', color='w', label='EXP Model Original Nodal Position', markerfacecolor='g', markersize=5),pll.Line2D([0], [0], marker='o', color='w', label='EXP Model Deformed Nodal Points', markerfacecolor='r', markersize=5),pll.Line2D([0], [0], marker='o', color='b', label='NUM Model Deformed Nodal Points', markerfacecolor='w', markersize=5)], loc = 'center')
        plt.subplot(221)
        plt.plot(orgx,orgy,'g.')
        plt.plot(x_new, y_new, 'r.')
        plt.plot(px,py, 'bo', fillstyle="none")
        plt.xlabel('X-Coordinate [mm]',fontsize=18)
        plt.ylabel('Y-Coordinate [mm]',fontsize=18)
        plt.grid("on")
        
        plt.subplot(222)
        plt.plot(x_new,dx_new, 'r.')
        plt.plot(x_new,x_rbf, 'bo', fillstyle="none")
        plt.xlabel('x [mm]',fontsize=18)
        plt.ylabel(r'$\Delta$x [mm]',fontsize=18)
        plt.grid("on")
        
        plt.subplot(223)
        plt.plot(y_new,dy_new, 'r.')
        plt.plot(y_new,y_rbf, 'bo', fillstyle="none")
        plt.xlabel('y [mm]',fontsize=18)
        plt.ylabel(r'$\Delta$y [mm]',fontsize=18)
        plt.grid("on")
        
        plt.subplot(224)
        plt.plot(z_new,dz_new, 'r.')
        plt.plot(z_new,z_rbf, 'bo', fillstyle="none")
        plt.xlabel('z [mm]',fontsize=18)
        plt.ylabel(r'$\Delta$z [mm]',fontsize=18)
        plt.grid("on")
        plt.subplots_adjust(left=0.08, right=0.94, bottom=0.08, top=0.90, wspace=0.2 ,hspace=0.55)
        fig6 = "Figure2_%s.pdf"%bb
        figpath6 = os.path.join(path1,fig6)
        plt.savefig(figpath6, dpi=600)
        plt.close(fig=fn)
        fn += 1

        # // Figure 3 plots the x, y and z displacements for the "NUM" model vs. "EXP" model
        plt.figure(fn,figsize=(12.0,8.5))
        plt.subplot(221)
        xx = np.linspace(min(dx_new),max(dx_new),100)
        yx = alpd[0]*xx + betd[0]
        plt.scatter(dx_new,x_rbf, s=15, c="b", alpha=0.5)
        plt.plot(xx,yx,'r')
        plt.title("X-Displacement")
        plt.xlabel(r'EXP Model $\Delta$x [mm]',fontsize=18)
        plt.ylabel(r'NUM Model $\Delta$x [mm]',fontsize=18)
        plt.grid("on")
        plt.text(min(dx_new), max(x_rbf)-0.2*(abs(max(x_rbf)-min(x_rbf))), "$R^2$ = %s, SEOE = %s\n RMS = %s"%(round(rsqid[0],4),round(seoedd[0],4),round(rmsd_f[0],4)),fontsize=14)
        
        plt.subplot(222)
        xy = np.linspace(min(dy_new),max(dy_new),100)
        yy = alpd[1]*xy + betd[1]
        plt.plot(xy,yy,'r')
        plt.scatter(dy_new,y_rbf,s=15, c="b", alpha=0.5)
        plt.title("Y-Displacement")
        plt.xlabel(r'EXP Model $\Delta$y [mm]',fontsize=18)
        plt.ylabel(r'NUM Model $\Delta$y [mm]',fontsize=18)
        plt.grid("on")
        plt.text(min(dy_new), max(y_rbf)-0.2*(abs(max(y_rbf)-min(y_rbf))), "$R^2$ = %s, SEOE = %s\n RMS = %s"%(round(rsqid[1],4),round(seoedd[1],4),round(rmsd_f[1],4)),fontsize=14)
        
        plt.subplot(223)
        xz = np.linspace(min(dz_new),max(dz_new),100)
        yz = alpd[2]*xz + betd[2]
        plt.plot(xz,yz,'r')
        plt.scatter(dz_new,z_rbf,s=15, c="b", alpha=0.5)
        plt.title("Z-Displacement")
        plt.xlabel(r'EXP Model $\Delta$z [mm]',fontsize=18)
        plt.ylabel(r'NUM Model $\Delta$z [mm]',fontsize=18)
        plt.grid("on")
        plt.text(min(dz_new), max(z_rbf)-0.2*(abs(max(z_rbf)-min(z_rbf))), "$R^2$ = %s, SEOE = %s\n RMS = %s"%(round(rsqid[2],4),round(seoedd[2],4),round(rmsd_f[2],4)),fontsize=14)
        plt.subplots_adjust(left=0.08, right=0.94, bottom=0.08, top=0.90, wspace=0.20 ,hspace=0.46)
        fig7 = "Figure3_%s.pdf"%bb
        figpath7 = os.path.join(path1,fig7)
        plt.savefig(figpath7, dpi=600)
        plt.close(fig=fn)
        fn += 1

        print('hey')

# ------------------------------------------------------------------------------------------------------------------
    # ---- Unindent the next lines of code in case two optimisation algorithms are used.
    # ---- Only unindent bedore the 2nd optimisation algorithm run, NOT during the 1st one
    
    # objls = sorted(objl, reverse=True)
    # xs = np.linspace(1,len(objls),len(objls))
    # objlsSP = sorted(objlSP, reverse=True)
    # xsSP = np.linspace(1,len(objlsSP),len(objlsSP))
    # plt.figure(1, figsize=(12.0,8.5))
    # plt.plot(xs,objls, "r*:")
    # plt.plot(xsSP,objlsSP, "bd:")
    # plt.xlabel('Optimisation Run',fontsize=18)
    # plt.ylabel('Optimisation Run Objective Function',fontsize=18)
    # plt.legend([algpat1,algpat2],fontsize=14)
    # plt.grid("on")
    # fig4 = "Figure4.pdf"
    # path3 = filepaths("path3")
    # figpath4 = os.path.join(path3,fig4)
    # plt.savefig(figpath4, dpi=600)
# ------------------------------------------------------------------------------------------------------------------

    return(best)

# b = post_process

