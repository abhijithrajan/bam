#! /usr/bin/env python 
import sys
import pyfits
from pyfits import getheader
from numpy import loadtxt, array
import numpy as np
from pylab import *
from os import listdir
import os
import time
import matplotlib
import math
import random
from scipy import stats

import unix as u
import param as par
import bin_calc as bc
import calc as c
import statistical as st

ion()

def get_airmass_JD(directory):
  global JD, imgs
  dir_contents = os.listdir(directory+'/'+par.data_sub_dir[0])
  imgs = [fn for fn in dir_contents if fn.startswith(par.image_name[0]) and fn.endswith(par.image_name[1])]
  imgs = sorted(imgs)
  headers = []
  airmass = []
  JD = []
  for i in range(len(imgs)):
    headers.append(pyfits.open(directory+'/'+par.data_sub_dir[0]+'/'+imgs[i]))
    JD.append(headers[i][0].header[par.JD_header[0]])
    airmass.append(headers[i][0].header[par.airmass[0]])
  return airmass, JD

def weightedmean(x,w):
	num=sum([x[i]*w[i] for i in range(len(w))])
	den=sum([w[i] for i in range(len(w))])
	return num/den

def convert(x):
	if 'INDEF' in x.strip(): x = -9999
	else: x = float(x.strip())
	return (x)

def plot(directory):
	global good_ref_stars,saturated_ref_stars,stat_reject_ref_stars
	good_ref_stars = []
	saturated_ref_stars = []
	stat_reject_ref_stars = []
	dissimilar_refs = []
	S_LC = []
	err_S_w_mean_norm = []
	sim_brightness_ref = []
	test_STD = []

	obj_name = directory.split('/')[-1]
	rootdir = par.rootdir[0]
	binned = par.binned_data[0]
	
	''' Finding the ds9 file '''
	dir_contents = os.listdir(directory+'/'+par.data_sub_dir[0])
	ds9s = [fn for fn in dir_contents if fn.startswith(par.ds9_name[0]) and fn.endswith(par.ds9_name[1])]
	ds9 = directory+'/'+par.data_sub_dir[0]+'/'+ds9s[0]
	ds9_file = loadtxt(ds9, unpack='true')
	xpos,ypos,fmax = ds9_file[0],ds9_file[1],ds9_file[4]

	''' Defining the number of reference stars '''
	fwhm = np.loadtxt(directory+'/'+par.data_sub_dir[0]+'/'+'fwhm.txt',usecols=(0,))
	N_ref_stars = len(ds9_file[0])-1
	N_ref_stars_orig = N_ref_stars
	
	bad_mag_refs = []
	mag_refs = [] 

	''' Read Julian Date and airmass info from .fits headers '''
	airmass, JD = get_airmass_JD(directory)
	airmass = np.array(airmass)
	JD = np.array(JD)
	turbulence_fwhm = np.ones(len(fwhm)) * np.array([random.random() for _ in xrange(len(fwhm))])*1e-12
	turbulence_airmass = np.ones(len(airmass)) * np.array([random.random() for _ in xrange(len(airmass))])*1e-12
	fwhm = fwhm+turbulence_fwhm		# A tiny pertubation is added as some of the FWHM values measured by IRAF are identical causing more FWHM elements to be revomved than what is desired.
	airmass = airmass+turbulence_airmass

	''' Calculate the seeing '''
	pixel_scale = par.pixel_scale[0]
	seeing = fwhm*pixel_scale

	''' Establishing directories '''
	dir_contents = os.listdir(directory)

	if par.variable_aperture[0]=="yes":
	  folders = [s for s in dir_contents if s.startswith("V") and s.endswith("M")]	  
	else:
	  folders = [s for s in dir_contents if s.startswith("A")]
	
	T = loadtxt(directory+'/'+folders[0]+'/xyf.list',converters={3: convert, 4:convert})

	''' Convert JD days --> hours '''
	JD_1 = JD.min()
	JD = (JD-JD_1)

	''' Creating 2D lists list[star][frame] '''
	flux_r = [[] for _ in range(N_ref_stars)]
	mag_r = [[] for _ in range(N_ref_stars)]
	err_r = [[] for _ in range(N_ref_stars)]
	err_r_MAG = [[] for _ in range(N_ref_stars)]
	flux_r_norm = [[] for _ in range(N_ref_stars)]
	flux_r_norm_all = [[] for _ in range(N_ref_stars)]
	err_r_norm = [[] for _ in range(N_ref_stars)]
	err_r_norm_all = [[] for _ in range(N_ref_stars)]
	stddev = [[] for _ in range(N_ref_stars)]
	stddev_norm = [[] for _ in range(N_ref_stars)]
	variance_norm = [[] for _ in range(N_ref_stars)]
	mag_refs = [[] for _ in range(N_ref_stars)]

	''' Obtaining the photometry data '''
	x, y, flux, mag, merr, msky = loadtxt(directory+'/'+folders[0]+'/xyf.list', unpack=True,converters={3: convert, 4:convert})
	x_target = x[::N_ref_stars+1]
	y_target = y[::N_ref_stars+1]


	''' Reading the Target Info '''
	flux_T = flux[::N_ref_stars+1]
	flux_T_MAG = mag[::N_ref_stars+1]
	
	''' Normalising Target Flux '''
	flux_target_norm = flux_T/np.median(flux_T)
	err_T = merr[::N_ref_stars+1]*flux_T/1.0857 # Convert from mag error to flux error
	err_T_MAG = merr[::N_ref_stars+1]+1e-6
	err_target_norm = err_T/np.median(flux_T)

	iteration = 0
	N_bad_refs = 0

	''' Reading reference star info and normalising '''
	for i in range(N_ref_stars):
	  flux_r[i] = np.array(flux[i+1::N_ref_stars+1])
	  mag_r[i] = np.array(mag[i+1::N_ref_stars+1])
	  err_r[i] = np.array(merr[i+1::N_ref_stars+1]*flux_r[i]/1.0857)
	  err_r_MAG[i] = np.array(merr[i+1::N_ref_stars+1] )
	  ''' Normalising the reference flux and errors '''
	  flux_r_norm[i] = np.array(flux_r[i]/np.median(flux_r[i]))
	  flux_r_norm_all[i] = np.array(flux_r[i]/np.median(flux_r[i]))
	  # Adding 1e-12 as IRAF has errors stopping at 0.0001 and I don't want to divide by zero:
	  err_r_norm[i] = np.array(err_r[i]/np.median(flux_r[i])+1e-12)
	  err_r_norm_all[i] = np.array(err_r[i]/np.median(flux_r[i])+1e-12)
	  stddev[i] = flux_r[i].std()
	  stddev_norm[i] = flux_r_norm[i].std()
	  variance_norm[i] = array(err_r_norm[i]*err_r_norm[i])
	  mag_refs[i] = array(np.median(mag_r[i]))

	while (iteration < 2):
	
	  REF_LCs = []
	  REF_LCs_median = []
	  err_ref_Rs = []
	  err_ref_Rs_median = []

	  delta_mag = []
	  
	  N_ref_stars = N_ref_stars-len(stat_reject_ref_stars)

	  print "\n\nITERATION: ",iteration
	  print 'Total Refs: ',N_ref_stars

	  ''' Removing bright/faint refs '''
	  if (iteration == 0):
	    for i in range(1,N_ref_stars+1): #Avoiding the 1st object in the list, which is the Target!
	      if (ds9_file[4][i] < par.cut_offs[4]):
	        print '\nRef : ',int(i),', Peak : ',ds9_file[4][i],u.yellow(', Too faint!')
	        saturated_ref_stars.append(i-1)    	
       
	      if (ds9_file[4][i] > par.cut_offs[5]):
	        print '\nRef: ',int(i),', Peak : ',ds9_file[4][i],u.red(', Too bright!')
	        saturated_ref_stars.append(i-1)

	    N_ref_stars = N_ref_stars - len(saturated_ref_stars)
	    print '\nRemoving bright/faint ref stars:',N_ref_stars,'refs remain'


	  ''' Selecting similar brightness refs  '''
	  if (iteration == 0 and N_ref_stars > par.cut_offs[2]):
	    diff = []
 	    sim_refs = []

 	    for i in range(N_ref_stars + len(saturated_ref_stars)):
 	      diff.append( [abs(np.median(mag_refs[i])-np.median(flux_T_MAG)),i] ) 	    
 	    diff = sorted(diff)

 	    if N_ref_stars > par.cut_offs[3]: # If the number is larger than the min number of ref stars
 	      for i in range(par.cut_offs[2]):
 	        if diff[i][1] not in saturated_ref_stars:
 	          sim_refs.append(diff[i][1])
 	    else:
 	      print 'Too many saturated or faint stars'
 	      sys.exit()
 	    
 	    dissimilar_refs = [x for x in np.arange(0,N_ref_stars_orig) if x not in sim_refs and x not in saturated_ref_stars]	    	
	    print 'Removing',len(dissimilar_refs),'refs which are too dissimilar'    

	  
	  ''' Creating lists with the good brightness refs '''
	  if (iteration == 0):
	    bad_ref_stars = saturated_ref_stars + dissimilar_refs
	    decent_refs = [x for x in np.arange(0,N_ref_stars_orig) if x not in bad_ref_stars]
 	    xpos,ypos,fmax = xpos[np.array(decent_refs)+1],ypos[np.array(decent_refs)+1],fmax[np.array(decent_refs)+1]

	    N_ref_stars = len(decent_refs)
	    print N_ref_stars,'Refs before statistical cut-offs'	
	    
	    flux_r = [flux_r[i] for i in decent_refs]
	    mag_r = [mag_r[i] for i in decent_refs]
	    err_r = [err_r[i] for i in decent_refs]
	    err_r_MAG = [err_r_MAG[i] for i in decent_refs]
	    flux_r_norm = [flux_r_norm[i] for i in decent_refs]
	    err_r_norm = [err_r_norm[i] for i in decent_refs]
	    stddev = [stddev[i] for i in decent_refs]
	    stddev_norm = [stddev_norm[i] for i in decent_refs]
	    variance_norm = [variance_norm[i] for i in decent_refs]
	    
	    mag_refs = [mag_refs[i] for i in decent_refs]
	  
	  ''' Removing the bad refs '''
	  if (iteration > 0):
	    flux_r = [flux_r[i] for i in good_ref_stars]
	    mag_r = [mag_r[i] for i in good_ref_stars]
	    err_r = [err_r[i] for i in good_ref_stars]
	    err_r_MAG = [err_r_MAG[i] for i in good_ref_stars]
	    flux_r_norm = [flux_r_norm[i] for i in good_ref_stars]
	    err_r_norm = [err_r_norm[i] for i in good_ref_stars]
	    stddev = [stddev[i] for i in good_ref_stars]
	    stddev_norm = [stddev_norm[i] for i in good_ref_stars]
	    variance_norm = [variance_norm[i] for i in good_ref_stars]
	    mag_refs = [mag_refs[i] for i in good_ref_stars]
	    xpos,ypos,fmax = xpos[good_ref_stars],ypos[good_ref_stars],fmax[good_ref_stars]
	    
	  ''' Find the mean of the raw flux of the reference stars '''
	  norm_ref_star_flux_sum = np.array([sum(a) for a in zip(*(flux_r_norm))])
	  norm_ref_star_flux_sum = norm_ref_star_flux_sum/N_ref_stars

	  ''' Initialising the 2D arrays to be used in the weighted mean calculations'''
	  weights = [[] for _ in range(N_ref_stars)]
	  weighted_mean = [[] for _ in range(N_ref_stars)]
	  norm_ref_star_flux_weight_sum_mean = [[] for _ in range(N_ref_stars)]

	  ''' Calculating the weighted mean '''
	  for i in range(N_ref_stars):
	    weights[i] = 1. / variance_norm[i]		# variance_norm = err_r_norm^2
	    weighted_mean[i] = weights[i]*flux_r_norm[i]	    
	  
	  weights_sum = np.array([mean(a) for a in zip(*(weights))])
	  weighted_mean_sum = np.array([mean(a) for a in zip(*(weighted_mean))])
	  weighted_mean_sum_norm = weighted_mean_sum/weighted_mean_sum.sum()	# Does sum to 1

	  ''' Creating Master Ref Curve with associated errors '''
	  norm_ref_star_flux_weight_sum_mean = weighted_mean_sum / weights_sum
	  sigma_weighted_mean_ref_stars = np.sqrt(1./(weights_sum))#*N_ref_stars)) #Erroneous error estimate

	  norm_ref_star_flux_median = np.array([np.median(a) for a in zip(*(flux_r_norm))])
	  sigma_median_ref_stars = np.sqrt(np.array([sum(a) for a in zip(*(variance_norm))])/(N_ref_stars-1)/N_ref_stars)

	  ''' Calculating the weighted mean excluding the ref star used
		  to create ref star light curves. i.e.
		  REF_LCs[2] = flux_r_norm[2]/norm_ref_star_flux_weight_sum_mean_Rs[2]
		  where norm_ref_star_flux_weight_sum_mean_Rs[2] does not include the
		  ref star with index 2.
		  '''	    
	  for i in range(N_ref_stars):
	    weights_Rs = weights
	    weighted_mean_Rs = weighted_mean
	    
	    weights_Rs = [x for x in weights_Rs if not (x == weights_Rs[i]).all()] # Removing the contribution by the comparison star
	    Rs = [x for x in flux_r_norm if not (x == flux_r_norm[i]).all()]
	    weights_sum_Rs = np.array([mean(a) for a in zip(*(weights_Rs))])
	    weighted_mean_Rs = [x for x in weighted_mean_Rs if not (x == weighted_mean_Rs[i]).all()]
	    
	    weighted_mean_sum_Rs = np.array([mean(a) for a in zip(*(weighted_mean_Rs))])	    
	    Rs_median = np.array([np.median(a) for a in zip(*(Rs))])
	      
	    weighted_mean_sum_norm_Rs = weighted_mean_sum_Rs/weighted_mean_sum_Rs.sum()
	    norm_ref_star_flux_weight_sum_mean_Rs = weighted_mean_sum_Rs / weights_sum_Rs

	  for i in range(N_ref_stars):	    
	    REF_LCs.append(flux_r_norm[i]/norm_ref_star_flux_weight_sum_mean_Rs)
	    REF_LCs_median.append(flux_r_norm[i]/Rs_median)

	    err_ref_Rs.append(np.sqrt(err_r_norm[i]**2 + np.median(sigma_weighted_mean_ref_stars)**2))	# IRAF errors
	    err_ref_Rs_median.append(np.sqrt(err_r_norm[i]**2 + np.median(sigma_median_ref_stars)**2))

	  master_ref = np.array([np.median(a) for a in zip(*(REF_LCs))])	  
	  master_ref = master_ref/np.median(master_ref)

	  master_ref_median = np.array([np.median(a) for a in zip(*(REF_LCs_median))])
	  master_ref_median = master_ref_median/np.median(master_ref_median)

	  ''' FINAL ERROR BARS using IRAF values '''
	  err_target_w_mean_norm = np.sqrt(err_target_norm**2 + np.median(sigma_weighted_mean_ref_stars)**2) # Target errors and error of combined ref LC added in quadrature.

	  err_target_median_norm = np.sqrt(err_target_norm**2 + np.median(sigma_median_ref_stars)**2) # Target errors and error of combined ref LC added in quadrature.


	  ''' Sum of reference star errors (as given by IRAF)  '''
	  err_r_norm_sum = np.array([np.median(a) for a in zip(*(err_r_norm))])/np.sqrt(N_ref_stars)		# Errors not using a weigthed mean.

	  ''' Check if in fact, the weighted mean has improved the STDDECV of LC. '''
	  if (iteration == 0):	
		if 	((flux_target_norm/norm_ref_star_flux_weight_sum_mean).std() <= (flux_target_norm/norm_ref_star_flux_sum).std()):
			print "\nThe weighted mean has ",u.green("improved")," the STDDEV of the LC.\n"
		else:
			print "Weighted mean ",u.red("not")," really helping here.\n"
			print (flux_target_norm/norm_ref_star_flux_weight_sum_mean).std(),"!<",(flux_target_norm/norm_ref_star_flux_sum).std(),"\n"
	
		if (sigma_weighted_mean_ref_stars.sum() <= err_r_norm_sum.sum() ):
			print "The weighted mean has ",u.green("improved")," the errors."
		else:
			print "Weighted mean ",u.red("not")," really helping the errors here."

	  
	  ''' The non detrended target light curve '''
	  T_LC = flux_target_norm/norm_ref_star_flux_weight_sum_mean
	  T_LC = T_LC/np.median(T_LC)

	  T_LC_med = flux_target_norm/norm_ref_star_flux_median
	  T_LC_med = T_LC_med/np.median(T_LC_med)
	
	  if (iteration == 0):
	    # Detrending AIRMASS
	    A,B,C = np.polyfit(airmass,master_ref, 2) # Fitting a 2D polynomial to the airmass
	    np_airmass = np.arange(airmass.min(),airmass.max(),1e-3)
	    poly2d_AIRMASS = A*np_airmass**2+B*np_airmass+C
	    airmass_effect =A*airmass**2+B*airmass+C# np.ones(len(airmass))#

	    A,B,C = np.polyfit(airmass,master_ref_median, 2) # Fitting a 2D polynomial to the airmass
	    np_airmass = np.arange(airmass.min(),airmass.max(),1e-3)
	    poly2d_AIRMASS = A*np_airmass**2+B*np_airmass+C
	    airmass_effect_median =A*airmass**2+B*airmass+C# np.ones(len(airmass))#
	 
	  if (iteration > 0 and par.variable_aperture[0]=="no"):   
	    # Detrending FWHM
	    slope_FWHM, intercept_FWHM, r_value_AIRMASS, p_value_AIRMASS, std_err_AIRMASS = stats.linregress(fwhm,master_ref)
	    print "AIRMASS params:\t",round(A,5),round(B,5),round(C,5)
	    print "FWHM params:\t",round(slope_FWHM,5),round(intercept_FWHM,5)
	    fwhm_effect = slope_FWHM*fwhm+intercept_FWHM

	    slope_FWHM, intercept_FWHM, r_value_AIRMASS, p_value_AIRMASS, std_err_AIRMASS = stats.linregress(fwhm,master_ref_median)
	    print "AIRMASS params:\t",round(A,5),round(B,5),round(C,5)
	    print "FWHM params:\t",round(slope_FWHM,5),round(intercept_FWHM,5)
	    fwhm_effect_median = slope_FWHM*fwhm+intercept_FWHM

	
	  ''' Creating the light curve for a test star.
	      The test star is amongst refs used for
	      the brown dwarf. '''	  
	  if (iteration > 0):
	    for k in range(N_ref_stars):
	      sim_brightness_ref.append(abs(np.median(mag_refs[k])-np.median(flux_T_MAG)))
	      Test_LC = REF_LCs[k]
	      test_STD.append([sim_brightness_ref[k],Test_LC.std(),k])
	    test_STD = sorted(test_STD)

	    n_1 = test_STD[0][2]
	    delta_mag = test_STD[0][0]
	    print "\nMost similar brightness ref: ",n_1+1," has Delta Mag = ",delta_mag,"with a STDDEV: ",round(test_STD[0][1],3)

	    T_LC = T_LC/airmass_effect
	    T_LC_med = T_LC_med/airmass_effect_median
	    S_LC = (flux_r_norm_all[n_1]/norm_ref_star_flux_weight_sum_mean_Rs)
	    S_LC = S_LC/airmass_effect
	    err_S_w_mean_norm = err_ref_Rs[n_1]
	    
	    if (par.variable_aperture[0]=="no"):
	      T_LC = T_LC/fwhm_effect
	      T_LC_med = T_LC_med/fwhm_effect_median
	      S_LC = Test_LC[n_1]/fwhm_effect
	    T_LC = T_LC/np.median(T_LC)
	    T_LC_med = T_LC_med/np.median(T_LC_med)
	    S_LC = S_LC/np.median(S_LC)#S_LC/np.mean(T_LC)
	    
	    for i in range(N_ref_stars):
	      REF_LCs[i] = REF_LCs[i]/airmass_effect
	      REF_LCs_median[i] = REF_LCs_median[i]/airmass_effect_median
	      if (par.variable_aperture[0]=="no"):
	        REF_LCs[i] = REF_LCs[i]/fwhm_effect
	        REF_LCs_median[i] = REF_LCs_median[i]/fwhm_effect_median
	      REF_LCs[i] = REF_LCs[i]/np.median(REF_LCs[i])
	      REF_LCs_median[i] = REF_LCs_median[i]/np.median(REF_LCs_median[i])

	    master_ref = np.array([np.median(a) for a in zip(*(REF_LCs))])
	    master_ref = master_ref/np.mean(master_ref)

	    master_ref_median = np.array([np.median(a) for a in zip(*(REF_LCs_median))])
	    master_ref_median = master_ref_median/np.median(master_ref_median)
	    
	    # Creating datafiles needed for plotting
	    if (iteration > 0 and binned == False):
	      f1 = open(rootdir+'/plots/data/unbinned/'+obj_name+'_unbinned_wm.txt', 'w+')
	      f2 = open(rootdir+'/plots/data/unbinned/'+obj_name+'_unbinned_med.txt', 'w+')
	      for i in range(len(T_LC)):
	        print >> f1, JD[i],T_LC[i],err_target_w_mean_norm[i],master_ref[i],err_r_norm_sum[i]
	        print >> f2, JD[i],T_LC_med[i],err_target_median_norm[i],master_ref_median[i],err_r_norm_sum[i]
	      f1.close()
	      f2.close()
	    else:
	      print rootdir+'/plots/data/binned/'+obj_name+'_wm.txt'
	      f1 = open(rootdir+'/plots/data/binned/'+obj_name+'_wm.txt', 'w+')
	      f2 = open(rootdir+'/plots/data/binned/'+obj_name+'_med.txt', 'w+')
	      for i in range(len(T_LC)):
	        print >> f1, JD[i],T_LC[i],err_target_w_mean_norm[i],master_ref[i],err_r_norm_sum[i],S_LC[i],err_S_w_mean_norm[i]
	        print >> f2, JD[i],T_LC_med[i],err_target_median_norm[i],master_ref_median[i],err_r_norm_sum[i]#,S_LC[i],err_S_w_mean_norm[i]
	      f1.close()	    
	      f2.close()	    
	  
	  c.calc(obj_name,N_ref_stars,N_bad_refs,airmass,np_airmass,poly2d_AIRMASS,\
fwhm,JD,T_LC,S_LC,err_target_w_mean_norm,err_S_w_mean_norm,norm_ref_star_flux_weight_sum_mean,\
REF_LCs,err_r_norm_sum,err_ref_Rs,master_ref,iteration,\
ds9_file,good_ref_stars,bad_ref_stars,stat_reject_ref_stars,N_ref_stars_orig,delta_mag,\
T_LC_med,err_target_median_norm,REF_LCs_median, err_ref_Rs_median, master_ref_median,\
xpos,ypos,fmax)
	  iteration += 1
