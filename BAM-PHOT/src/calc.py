import numpy as np
from itertools import izip_longest
from pylab import *
import math
from scipy import stats

import unix as u
import statistical as st
import param as par

import datetime
import os

def plot_them(iteration,obj_name,JD,T_LC,S_LC,err_target_w_mean_norm,err_S_w_mean_norm,\
REF_LCs,master_ref,err_r_norm_sum,airmass,np_airmass,poly2d_AIRMASS,N_ref_stars,\
err_ref_Rs,fwhm,JD_unbinned,T_LC_unbinned,err_target_w_mean_norm_unbinned,\
master_ref_unbinned,err_r_norm_sum_unbinned,JD_binned,T_LC_binned,\
T_LC_UNCERTAINTY,T_LC_binned_uncertainty,T_LC_binned_IRAF_uncertainty,dividers,binned,\
T_LC_med,err_target_median_norm,REF_LCs_median, err_ref_Rs_median, master_ref_median):
 
  rainbow_colors = iter(cm.rainbow(np.linspace(0, 1, N_ref_stars)))
  
  if (iteration == 0):
    target_colour = '#D8D8D8'
    master_ref_colour = '#D8D8D8'
  else:
    target_colour = 'red'
    master_ref_colour = 'red'
 
  #'''
  # FIGURE 1:
  #---------------------------------------------------------------------------------------
  fig=plt.figure(1,figsize=(8.27,11.69))
  
  # Light curve plot
  main_panel=fig.add_subplot(3,2,1)
  main_panel.cla()
#  if (iteration > 0):
#    plt.errorbar(JD,S_LC,yerr=err_S_w_mean_norm,fmt='.',elinewidth=1.3,\
#    capsize=3,markersize=12,markeredgecolor='black',color="blue", ecolor="black")
  plt.errorbar(JD*24.,T_LC,yerr=err_target_w_mean_norm,fmt='.',elinewidth=1.3,\
  capsize=3,markersize=12,markeredgecolor='black',color=target_colour, ecolor="black")
  ylabel('Relative Flux')

  # Master ref plot
  middle_panel=fig.add_subplot(3,2,3,sharex=main_panel, sharey=main_panel)
  middle_panel.cla()
  for i in range(len(REF_LCs)):
    plt.plot(JD*24.,REF_LCs[i],'.k',alpha=0.1)
  plt.errorbar(JD*24.,master_ref,yerr=err_r_norm_sum,fmt='.',elinewidth=1.3,\
  capsize=3,markersize=12,markeredgecolor='black',color=master_ref_colour, ecolor="black")
  ylabel('Relative Flux')
  #plt.ylim(T_LC_binned.min()-0.01,T_LC_binned.max()+0.01)
  
  # Airmass trend plot
  if (iteration == 0):
    bottom_panel=fig.add_subplot(3,2,5)
    bottom_panel.cla()
    plt.plot(airmass,master_ref,'ok')
    plt.plot(np_airmass,poly2d_AIRMASS,'-r')
    xlabel('Airmass')
    ylabel('Relative Flux')
  else:
    bottom_panel=fig.add_subplot(3,2,5)
    bottom_panel.cla()
    plt.plot(JD*24.,fwhm,'ok')
    xlabel('Time (Hours)')
    ylabel('FWHM')
    
  # Ref star plots on the right hand side
  
  # Dirty code:
  #os.system('rm '+par.rootdir[0]+'/plots/data/binned/'+obj_name+'_refs.txt')
  
  if (iteration == 0):
    for i in range(N_ref_stars):
      side_panel = fig.add_subplot(N_ref_stars,2,2*(i+1),sharex=main_panel, sharey=main_panel)
      plt.errorbar(JD*24., REF_LCs[i],yerr=err_ref_Rs[i],fmt='.',elinewidth=1.3,capsize=3,markersize=12,ecolor='black',color=next(rainbow_colors),markeredgecolor='black')

  else:
    for i in range(N_ref_stars):
      side_panel = fig.add_subplot(N_ref_stars,2,2*(i+1),sharex=main_panel, sharey=main_panel)
      side_panel.cla() 
      plt.errorbar(JD*24., REF_LCs[i],yerr=err_ref_Rs[i],fmt='.',elinewidth=1.3,capsize=3,markersize=12,ecolor='black',color=next(rainbow_colors),markeredgecolor='black')
      #plt.ylim(np.array([T_LC_binned.min(),master_ref.min()]).min()-0.01,np.array([T_LC_binned.max(),master_ref.max()]).max()+0.01)
         
      
      if binned == False: f = open(par.rootdir[0]+'/plots/data/unbinned/'+obj_name+'_refs_wm.txt', 'a+')
      if binned == True: f = open(par.rootdir[0]+'/plots/data/binned/'+obj_name+'_refs_wm.txt', 'a+')
      print >> f,'----------- REF ',i+1,' -----------'
      for j in range(len(REF_LCs[i])):
        print >> f, JD[j], REF_LCs[i][j],err_ref_Rs[i][j]
      print >> f,'\n'
      f.close() 
      
  
  if par.variable_aperture[0] == "yes":
    aperture = "variable"
  else:
    aperture = "constant"
  minorticks_on()
#  plt.ylim(0.90,1.10)
  plt.savefig(par.rootdir[0]+'/plots/'+str(iteration)+'_'+obj_name+'_'+aperture+'_aperture.pdf')
  #plt.draw()
  if (iteration == 0):
    clf()
  #---------------------------------------------------------------------------------------


  #'''
  # FIGURE 3:
  #---------------------------------------------------------------------------------------
  fig3=plt.figure(3,figsize=(8.27,11.69))
  
  # Light curve plot
  main_panel=fig3.add_subplot(3,2,1)
  main_panel.cla()
#  if (iteration > 0):
#    plt.errorbar(JD,S_LC,yerr=err_S_w_mean_norm,fmt='.',elinewidth=1.3,\
#    capsize=3,markersize=12,markeredgecolor='black',color="blue", ecolor="black")
  plt.errorbar(JD*24.,T_LC_med,yerr=err_target_median_norm,fmt='.',elinewidth=1.3,\
  capsize=3,markersize=12,markeredgecolor='black',color=target_colour, ecolor="black")
  ylabel('Relative Flux')

  # Master ref plot
  middle_panel=fig3.add_subplot(3,2,3,sharex=main_panel, sharey=main_panel)
  middle_panel.cla()
  for i in range(len(REF_LCs_median)):
    plt.plot(JD*24.,REF_LCs_median[i],'.k',alpha=0.1)
  plt.errorbar(JD*24.,master_ref_median,yerr=err_r_norm_sum,fmt='.',elinewidth=1.3,\
  capsize=3,markersize=12,markeredgecolor='black',color=master_ref_colour, ecolor="black")
  ylabel('Relative Flux')
  #plt.ylim(T_LC_binned.min()-0.01,T_LC_binned.max()+0.01)
  
  # Airmass trend plot
  if (iteration == 0):
    bottom_panel=fig3.add_subplot(3,2,5)
    bottom_panel.cla()
    plt.plot(airmass,master_ref_median,'ok')
    plt.plot(np_airmass,poly2d_AIRMASS,'-r')
    xlabel('Airmass')
    ylabel('Relative Flux')
  else:
    bottom_panel=fig3.add_subplot(3,2,5)
    bottom_panel.cla()
    plt.plot(JD*24.,fwhm,'ok')
    xlabel('Time (Hours)')
    ylabel('FWHM')
    
  # Ref star plots on the right hand side
  rainbow_colors = iter(cm.rainbow(np.linspace(0, 1, N_ref_stars)))
  if (iteration == 0):
    for i in range(N_ref_stars):
      side_panel = fig3.add_subplot(N_ref_stars,2,2*(i+1),sharex=main_panel, sharey=main_panel)
      plt.errorbar(JD*24., REF_LCs_median[i],yerr=err_ref_Rs_median[i],fmt='.',elinewidth=1.3,capsize=3,markersize=12,ecolor='black',color=next(rainbow_colors),markeredgecolor='black')

  else:
    for i in range(N_ref_stars):
      side_panel = fig3.add_subplot(N_ref_stars,2,2*(i+1),sharex=main_panel, sharey=main_panel)
      side_panel.cla() 
      plt.errorbar(JD*24., REF_LCs_median[i],yerr=err_ref_Rs_median[i],fmt='.',elinewidth=1.3,capsize=3,markersize=12,ecolor='black',color=next(rainbow_colors),markeredgecolor='black')
      #plt.ylim(np.array([T_LC_binned.min(),master_ref.min()]).min()-0.01,np.array([T_LC_binned.max(),master_ref.max()]).max()+0.01)
         
      
      if binned == False: f = open(par.rootdir[0]+'/plots/data/unbinned/'+obj_name+'_refs_med.txt', 'a+')
      if binned == True: f = open(par.rootdir[0]+'/plots/data/binned/'+obj_name+'_refs_med.txt', 'a+')
      print >> f,'----------- REF ',i+1,' -----------'
      for j in range(len(REF_LCs[i])):
        print >> f, JD[j], REF_LCs[i][j],err_ref_Rs[i][j]
      print >> f,'\n'
      f.close() 
  

  if par.variable_aperture[0] == "yes":
    aperture = "variable"
  else:
    aperture = "constant"
  minorticks_on()
#  plt.ylim(0.90,1.10)
  plt.savefig(par.rootdir[0]+'/plots/'+str(iteration)+'_'+obj_name+'_'+aperture+'_aperture_median.pdf')
    #plt.draw()
  if (iteration == 0):
    clf()
  #---------------------------------------------------------------------------------------


  # FIGURE 2:
  #---------------------------------------------------------------------------------------
  if binned == False:
    fig2=plt.figure(2,figsize=(11.69,8.27))
    top_panel=fig2.add_subplot(2,1,1)
    top_panel.cla()
    plt.plot(JD_unbinned, master_ref_unbinned,'.',color='#FF9700',markersize=12,markeredgecolor='black')
    plt.plot(JD_unbinned,T_LC_unbinned,'.',color='#0C5AA6',alpha=1.0,markersize=16,markeredgecolor='black') # Unbinned data
    #plt.errorbar(JD_binned,T_LC_binned,yerr=T_LC_UNCERTAINTY,fmt='.',elinewidth=1.3,\
    #capsize=3,markersize=12,markeredgecolor='black',color='red', ecolor="black",alpha=1.0)

    for i in range(len(dividers)):
      top_panel.axvline(dividers[i],linewidth=2, ls='--',color='r',alpha=0.1)							# Divider lines
    xlabel('Time (Hours)')
    ylabel('Relative Flux')
    #plt.ylim(T_LC_binned.min()-0.05,T_LC_binned.max()+0.05)
    #plt.ylim(0.95,1.06)
    minorticks_on()
    plt.savefig(par.rootdir[0]+'/plots/'+obj_name+'_LC2.pdf')
    plt.show()
    #---------------------------------------------------------------------------------------

def calc(obj_name,N_ref_stars,N_bad_refs,airmass,np_airmass,poly2d_AIRMASS,\
fwhm,JD,T_LC,S_LC,err_target_w_mean_norm,err_S_w_mean_norm,norm_ref_star_flux_weight_sum_mean,\
REF_LCs,err_r_norm_sum,err_ref_Rs,master_ref,iteration,\
ds9_file,good_ref_stars,bad_ref_stars,stat_reject_ref_stars,N_ref_stars_orig,delta_mag,\
T_LC_med,err_target_median_norm,REF_LCs_median, err_ref_Rs_median, master_ref_median,\
xpos,ypos,fmax):
   
  ''' Photometric quality and amplitude calculation '''
  phot_qual = []
  for i in range(len(T_LC)):
    phot_qual.append(np.sqrt(err_target_w_mean_norm[i]**2+err_r_norm_sum[i]**2))
    if T_LC[i] == T_LC.max():
      top_err = err_target_w_mean_norm[i]
    if T_LC[i] == T_LC.min():
      bottom_err = err_target_w_mean_norm[i]
  amplitude = (T_LC.max()-T_LC.min())*100.
  amplitude_err = np.sqrt(top_err**2+bottom_err**2)*100.
  
  binned = par.binned_data[0]	# See if the param file says if images have been combined
  JD_unbinned = []
  T_LC_unbinned = []
  err_target_w_mean_norm_unbinned = []
  master_ref_unbinned = []
  err_r_norm_sum_unbinned = []
  JD_binned = []
  T_LC_binned = []
  T_LC_UNCERTAINTY = []
  dividers = []
  T_LC_binned_uncertainty = []
  T_LC_binned_IRAF_uncertainty = []
  
  if iteration > 0 and binned == False:
    JD_unbinned,T_LC_unbinned,err_target_w_mean_norm_unbinned,master_ref_unbinned,err_r_norm_sum_unbinned = loadtxt(par.rootdir[0]+'/plots/data/unbinned/'+obj_name+'_unbinned_wm.txt', unpack=True)
  
    # Binning Data
    divide = []
    index_divide = [0]
    for i in range(len(JD_unbinned)-2):
	  if (JD_unbinned[i+2]-JD_unbinned[i+1]) >= (JD_unbinned[i+1]-JD_unbinned[i])*4.0:
	    index_divide.append(i+2)
	    divide.append(JD_unbinned[i+1] + (JD_unbinned[i+2]-JD_unbinned[i+1])/2.)

    index_divide.append(i+2)
    divide_at = (JD_unbinned[-1] - JD_unbinned[0])/5.
    time_step = np.arange(JD_unbinned[0],JD_unbinned[-1],divide_at)
    k = 1
    if len(index_divide) < 3:
      index_divide = [0]
      for i in range(len(JD_unbinned-2)):
        if (divide_at*k <= JD_unbinned[i] <= divide_at*(k+1)) and (k < len(time_step)-1):
          #print divide_at*k,JD_unbinned[i],divide_at*(k+1)
          index_divide.append(i)
          divide.append(time_step[k])
          k += 1
		
	# SPECIAL CASES 
	if obj_name[3:5] == '07':
	  del index_divide[4]
	  index_divide.append(80)
	  index_divide.append(100)
	  index_divide.append(120)
	  index_divide.append(140)
	  index_divide.append(160)
	  index_divide.append(180)
	  index_divide.append(200)
	  index_divide.append(217)
	  divide.append(JD_unbinned[80+1] + (JD_unbinned[80+2]-JD_unbinned[80+1])/2.)
	  divide.append(JD_unbinned[100+1] + (JD_unbinned[100+2]-JD_unbinned[100+1])/2.)
	  divide.append(JD_unbinned[120+1] + (JD_unbinned[120+2]-JD_unbinned[120+1])/2.)
	  divide.append(JD_unbinned[140+1] + (JD_unbinned[140+2]-JD_unbinned[140+1])/2.)
	  divide.append(JD_unbinned[160+1] + (JD_unbinned[160+2]-JD_unbinned[160+1])/2.)
	  divide.append(JD_unbinned[180+1] + (JD_unbinned[180+2]-JD_unbinned[180+1])/2.)
	  divide.append(JD_unbinned[200+1] + (JD_unbinned[200+2]-JD_unbinned[200+1])/2.)
	  divide.append(JD_unbinned[217+1] + (JD_unbinned[217+2]-JD_unbinned[217+1])/2.)

	if obj_name[3:5] == '12':
	  del index_divide[-1]
	  del divide[-1]

	if obj_name[3:5] == '17':
	  del index_divide[-1]
	  del divide[-1]

	if obj_name[3:5] == '37':
	  del index_divide[-1]
	  del divide[-1]

	if obj_name[3:5] == '46':
	  del index_divide[5]
	  del divide[5]

    JD_bin = []
    JD_binned = []
    T_LC_bin = []
    T_LC_binned = []
    
    err_target_w_mean_norm_bin = []
    master_bin = []
    master_binned = []
    master_binned_uncertainty = []
    master_binned_IRAF_uncertainty = []
    err_r_norm_sum_bin = []
 
    for j in range (len(index_divide)-1): # Since we start counting at 0	
    	dividers.append(divide[j-1]) # Used to draw the vertical lines
	
  	JD_bin.append(np.array(JD_unbinned[index_divide[j]:index_divide[j+1]]))	# Define which points are within the bin
  	JD_binned.append(np.median(JD_bin[j]))					# Calculate mean of points within the bin
  	
  	T_LC_bin.append(T_LC_unbinned[index_divide[j]:index_divide[j+1]])
  	err_target_w_mean_norm_bin.append(err_target_w_mean_norm_unbinned[index_divide[j]:index_divide[j+1]])
  	T_LC_binned.append(np.median(T_LC_bin[j]))
  	T_LC_binned_uncertainty.append(st.sigma_clip(T_LC_bin[j], 3.0, None, maout=True).std()/np.sqrt(len(T_LC_bin[j])))
  	T_LC_binned_IRAF_uncertainty.append(err_target_w_mean_norm_bin[j].mean()/np.sqrt(len(err_target_w_mean_norm_bin[j])))
  		
  	master_bin.append(master_ref_unbinned[index_divide[j]:index_divide[j+1]])
  	err_r_norm_sum_bin.append(master_ref_unbinned[index_divide[j]:index_divide[j+1]])
  	master_binned.append(np.median(master_bin[j]))
  	master_binned_uncertainty.append(st.sigma_clip(master_bin[j], 3.0, None, maout=True).std()/np.sqrt(len(master_bin[j])))
  	master_binned_IRAF_uncertainty.append(err_r_norm_sum_bin[j].mean()/np.sqrt(len(err_r_norm_sum_bin[j])))
  
    T_LC_UNCERTAINTY = []
    for i in range(len(T_LC_binned)):
      if T_LC_binned_uncertainty[i] > T_LC_binned_IRAF_uncertainty[i]:
        T_LC_UNCERTAINTY.append(T_LC_binned_uncertainty[i])
      else:
        T_LC_UNCERTAINTY.append(T_LC_binned_IRAF_uncertainty[i])
  
    master_UNCERTAINTY = []
    for i in range(len(T_LC_binned)):
      if master_binned_uncertainty[i] > master_binned_IRAF_uncertainty[i]:
        master_UNCERTAINTY.append(master_binned_uncertainty[i])
      else:
        master_UNCERTAINTY.append(master_binned_IRAF_uncertainty[i])
      
    T_LC_binned = np.array(T_LC_binned)
    T_LC_binned_uncertainty = np.array(T_LC_binned_uncertainty)
    T_LC_binned_IRAF_uncertainty = np.array(T_LC_binned_IRAF_uncertainty)
    T_LC_UNCERTAINTY = np.array(T_LC_UNCERTAINTY)
  
    master_binned_uncertainty = np.array(master_binned_uncertainty)
    master_UNCERTAINTY = np.array(master_UNCERTAINTY)
  
  MAD_factor = par.cut_offs[1]
  STD_factor = par.std_cutoff[0]
  MAD_cutoff = []
  STD_cutoff = []
   
  for i in range(len(REF_LCs)):
    MAD_cutoff.append(st.MAD(REF_LCs[i]))
    STD_cutoff.append(REF_LCs[i].std())
  
  MAD_median = np.median(MAD_cutoff)
  STD_median = np.mean(STD_cutoff) # changed from median to mean
  MAD_STD = np.std(MAD_cutoff)
  STD_STD = np.std(STD_cutoff)
  
  ref_std = []
  for i in range(len(REF_LCs)):
    ref_std.append(REF_LCs[i].std())
  ref_std = np.array(ref_std)
  mean_ref_std = ref_std.mean()

  robust = st.robust(T_LC,err_target_w_mean_norm)
  chi2m = st.chi2_master_ref(T_LC,err_target_w_mean_norm,master_ref,err_r_norm_sum)
  #robust_binned = st.robust(T_LC_binned,T_LC_binned_uncertainty)
  #chi2_binned = st.chi2_master_ref(T_LC_binned,T_LC_UNCERTAINTY,master_binned,master_UNCERTAINTY)
  
  # Print to screen
  print "\n"
  print "\tChi2m:\t\t\t",round(chi2m,2)#,round(chi2_binned,2)
  print "\tRobust:\t\t\t",round(robust,2)#,round(robust_binned,2)
  #print "Ave. Uncertainty Target (binned): ",T_LC_UNCERTAINTY.mean()
  print "\tMean Ref std:\t\t",round(mean_ref_std,5)
  print "\tSTDDEV Target:\t\t",round(T_LC.std(),5)#,round(T_LC_binned.std(),5)
  print "\tSTDDEV Master Ref:\t",round(master_ref.std(),5)
  print "\tAmplitude:\t\t",round(amplitude,2),"+-",round(amplitude_err,2),"\n"
   
  if (iteration == 0):
    print "\tRef Stars:\t\t","cut-off: "+str(round(STD_median+STD_factor*STD_STD,5))+"\n"
    print "\tChi2:","\t\t","Robust:","\t","STD:","\t\t","xy-pos:","\t","Peak Flux:","\t","Ref Number:"
    print "\t--------------------------------------------------------------------------------------------------------"
        
    for i in range(N_ref_stars):
      if REF_LCs[i].std() > STD_median+STD_factor*STD_STD and N_ref_stars > par.cut_offs[3]:
        print "\t",round(st.chi2(REF_LCs[i],err_ref_Rs[i]),2),\
        '\t\t',round(st.robust(REF_LCs[i],err_ref_Rs[i]),2),\
        '\t\t',u.red(str(round(REF_LCs[i].std(),5))),\
        '\t',str(xpos[i]),',',str(ypos[i]),\
        '\t',str(fmax[i])
        stat_reject_ref_stars.append(i)

      elif st.robust(REF_LCs[i],err_ref_Rs[i]) >= 100. and N_ref_stars > par.cut_offs[3]:
        print "\t",round(st.chi2(REF_LCs[i],err_ref_Rs[i]),2),\
        '\t\t',u.red(str(round(st.robust(REF_LCs[i],err_ref_Rs[i]),2))),\
        '\t\t',round(REF_LCs[i].std(),5),\
        '\t',str(xpos[i]),',',str(ypos[i]),\
        '\t',str(fmax[i])
        stat_reject_ref_stars.append(i)

      elif st.chi2(REF_LCs[i],err_ref_Rs[i]) >= 50. and N_ref_stars > par.cut_offs[3]:
        print "\t",round(st.chi2(REF_LCs[i],err_ref_Rs[i]),2),\
        '\t\t',u.red(str(round(st.robust(REF_LCs[i],err_ref_Rs[i]),2))),\
        '\t\t',round(REF_LCs[i].std(),5),\
        '\t',str(xpos[i]),',',str(ypos[i]),\
        '\t',str(fmax[i])
        stat_reject_ref_stars.append(i)
       
      else:
      	good_ref_stars.append(i)

        print "\t",round(st.chi2(REF_LCs[i],err_ref_Rs[i]),2),\
        '\t\t',round(st.robust(REF_LCs[i],err_ref_Rs[i]),2),\
        '\t\t',round(REF_LCs[i].std(),5),u.green('OK'),\
        '\t',str(xpos[i]),',',str(ypos[i]),\
        '\t',str(fmax[i]),'\t\t',int(i+1.)
    print "\t--------------------------------------------------------------------------------------------------------" 

  if (iteration > 0):

    chi2_S = st.chi2_master_ref(S_LC,err_S_w_mean_norm,master_ref,err_r_norm_sum)
    
    DOF = len(master_ref)-1
    chance_variable = round((1-stats.chi2.cdf(chi2m*DOF, DOF))*100.,4)
    chance_variable_S = round((1-stats.chi2.cdf(chi2_S*DOF, DOF))*100.,4)

    print "\tp-val:\t\t\t",str(round(chance_variable,2))+"\n"
    print "\tRef Stars:\t\t","cut-off: "+str(round(STD_median+STD_factor*STD_STD,5))+"\n"
    print "\tChi2:","\t\t","Robust:","\t","STD:","\t\t","p-val:","\t\t","xy-pos:","\t\t","Peak Flux:"
    print "\t--------------------------------------------------------------------------------------------------------"
        
    for i in range(len(REF_LCs)):
      DOF = len(master_ref)-1
      chi2ref = st.chi2_master_ref(REF_LCs[i],err_ref_Rs[i],master_ref,err_r_norm_sum)
      print "\t",round(st.chi2(REF_LCs[i],err_ref_Rs[i]),2),\
        '\t\t',round(st.robust(REF_LCs[i],err_ref_Rs[i]),2),\
        '\t\t',u.blue(str(round(REF_LCs[i].std(),5))),\
        '\t',u.green(str(round((1-stats.chi2.cdf(chi2ref*DOF, DOF))*100.,4))),\
        '\t\t',str(xpos[i]),',',str(ypos[i]),\
        '\t',str(fmax[i])
    print "\t--------------------------------------------------------------------------------------------------------"
  
    if (T_LC.std() > master_ref.std() and robust >= 1.0 and chance_variable <= 5.0):
      variable = 'y'    
    else:
      variable = 'n'

    f = open(par.rootdir[0]+'/data.txt', 'a+')
    print >> f, obj_name,\
    '\t',round(chi2m,2),\
    '\t',round(robust,2),\
    '\t',round(chance_variable,2),\
    '\t',round(np.median(phot_qual),5),\
    '\t',round(mean_ref_std,6),\
    '\t',round(master_ref.std(),6),\
    '\t',str(len(good_ref_stars)+len(stat_reject_ref_stars)+len(bad_ref_stars)),\
    '\t',str(len(stat_reject_ref_stars)),\
    '\t',str(len(good_ref_stars)),\
    '\t',DOF,\
    '\t',round(amplitude,2),\
    '\t',round(amplitude_err,2),\
    '\t',round((JD[-1]-JD[0])*24.,2),\
    '\t',chance_variable,\
    '\t',variable,\
    '\t',chance_variable_S,\
    '\t',N_ref_stars_orig,\
    '\t',delta_mag,\
    '\t',T_LC.std()
    f.close()

    if binned == False: f = open(par.rootdir[0]+'/plots/data/unbinned/%s_sim_ref.txt'%(obj_name), 'a+')
    if binned == True: f = open(par.rootdir[0]+'/plots/data/binned/%s_sim_ref.txt'%(obj_name), 'a+')
    for j in range(len(S_LC)):
      print >> f, JD[j], S_LC[j],err_S_w_mean_norm[j]
    f.close()	   

  
  #'''
  plot_them(iteration,obj_name,JD,T_LC,S_LC,err_target_w_mean_norm,\
  err_S_w_mean_norm,REF_LCs,master_ref,err_r_norm_sum,airmass,np_airmass,poly2d_AIRMASS,\
  N_ref_stars,err_ref_Rs,fwhm,JD_unbinned,T_LC_unbinned,err_target_w_mean_norm_unbinned,\
  master_ref_unbinned,err_r_norm_sum_unbinned,JD_binned,T_LC_binned,T_LC_UNCERTAINTY,T_LC_binned_uncertainty,T_LC_binned_IRAF_uncertainty,dividers,binned,\
  T_LC_med,err_target_median_norm,REF_LCs_median, err_ref_Rs_median, master_ref_median)
  #'''
