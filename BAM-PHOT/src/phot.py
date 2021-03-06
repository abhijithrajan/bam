#! /usr/bin/env python 

#=======================================#
#  BAM-Phot
#  Originally written to work with the NTT.
#  Updated to work with a range of instruments
#  
#  Written by: Paul Anthony Wilson
#  paw@astro.ex.ac.uk
#  Updated by: Abhijith Rajan
#  arajan6@asu.edu
#=======================================#

import os, sys
import shutil
from os import getcwd,listdir,chdir
import numpy as np
from pylab import loadtxt

import pyfits
from pyraf import iraf

import param as par
import unix as u

def load_phot():
   iraf.noao(_doprint=0)
   iraf.digiphot(_doprint=0)
   iraf.daophot(_doprint=0)

def set_phot_params(fwhm,newdir):
   # Set centerpars
   iraf.centerpars.setParam('calgori','centroid')
   iraf.centerpars.setParam('cbox',5)
   iraf.centerpars.setParam('cmaxite',10)
   iraf.centerpars.setParam('clean','yes')
   iraf.centerpars.setParam('maxshift',5)

   # Set datapars
   iraf.datapars.setParam('fwhmpsf',fwhm)
   iraf.datapars.setParam('sigma',7.8)
   #iraf.datapars.setParam('datamin',1000)
   #iraf.datapars.setParam('datamax',45000)
   iraf.datapars.setParam('epadu',par.gain[0])
   iraf.datapars.setParam('readnoi',par.read_noise[0])
   iraf.datapars.setParam('exposure','EXPTIME')
   #iraf.datapars.setParam('airmass','HIERARCH ESO TEL AIRM START')
   #iraf.datapars.setParam('filter','HIERARCH ESO INS FILT1 NAME')
   #iraf.datapars.setParam('obstime','READTIME')

   # Set findpars
   iraf.findpars.setParam('thresho',3) # Threshold in sigma for feature detection
   iraf.findpars.setParam('nsigma',1.5) # Width of convolution kernel in sigma - whatever that means...

   # Set fitskypars
   iraf.fitskypars.setParam('salgorithm','mode')
   #iraf.fitskypars.setParam('annulus',12) # Inner radius of sky annulus in scale units
   #iraf.fitskypars.setParam('dannulus',15) # Width of sky annulus in scale units

   # Saving parameters in files
   #iraf.centerpars.saveParList(filename=newdir+'/params/center.par')
   #iraf.datapars.saveParList(filename=newdir+'/params/data.par')
   #iraf.fitskypars.saveParList(filename=newdir+'/params/fitsky.par')
   #iraf.photpars.saveParList(filename=newdir+'/params/phot.par')
   
   # Saving parameters in files
   iraf.centerpars.saveParList(filename='center.par')
   iraf.datapars.saveParList(filename='data.par')
   iraf.fitskypars.saveParList(filename='fitsky.par')
   iraf.photpars.saveParList(filename='phot.par')

   # Set some IRAF settings
   iraf.set(imtype="fits")


def phot(directory,folders):

  load_phot()
  
  print '____________________________________________________'

  dir_contents = os.listdir(directory+'/'+par.data_sub_dir[0])
  ds9s = [fn for fn in dir_contents if fn.startswith('ds9') and fn.endswith('v2.reg')]
  ds9 = directory+'/'+par.data_sub_dir[0]+'/'+ds9s[0]

  imgs = [fn for fn in dir_contents if fn.startswith(par.image_name[0]) and fn.endswith(par.image_name[1])]
  imgs = sorted(imgs)

  # Determine median FWHM for this dataset
  fwhm = np.loadtxt(directory+'/'+par.data_sub_dir[0]+'/fwhm.txt',usecols=(0,))

  median_fhwm = np.median(fwhm)
  seeing_value = median_fhwm*par.pixel_scale[0]
  print "\nMedian seeing: "+str(round(seeing_value,2))+"\""
  
  variable_aperture=par.variable_aperture[0]

  if variable_aperture=="no":
    Ap1 = 1.0*median_fhwm
    Ap2 = 1.0*median_fhwm
    Ap3 = 1.0*median_fhwm 
    Ap4 = 1.0*median_fhwm
  else:
    Ap1 = 1.0
    Ap2 = 1.0
    Ap3 = 2.0
    Ap4 = 2.0
  
  Ap5 = 10
  Ap6 = 10

  if Ap1 == Ap2:
	a_range = [Ap1]
	skyrad = [Ap3]
  else:
	a_range = np.arange(float(Ap1),float(Ap2)+0.25,0.25) 	# Step size
	skyrad=np.arange(float(Ap3),float(Ap4)+1.0,1.0)
  anwidth=np.arange(float(Ap5),float(Ap6)+5.0,5.0)

  for j in range(0,len(a_range)):
	 for k in range(0,len(skyrad)):
		 for l in range(0,len(anwidth)):
			 if variable_aperture=="no":
			   newsubdir = "/A" + str(a_range[j])
			   os.system('rm -r '+directory+'/A*')
			 else:
			   newsubdir = "/V" + str(a_range[j]) + "xFWHM"
			   #os.system('rm -r '+directory+'/V*')
			 
			 newdir=directory+newsubdir
			 
			 if os.path.isdir(newdir):
			   shutil.rmtree(newdir)
			 u.mkdir(newdir)
			 
			 mag_dir = newdir+"/mag"
			 params = newdir+"/params"
			 
			 u.mkdir(mag_dir)
			 u.mkdir(params)    
			 	
			 # Set outer photpars
			 iraf.fitskypars.setParam('annulus',skyrad[k]) # Inner radius of sky annulus in scale units
			 iraf.fitskypars.setParam('dannulus',anwidth[l]) # Width of sky annulus in scale units

			 iraf.phot.setParam('interactive','no')
			 iraf.phot.setParam('verify','no')
			 iraf.phot.setParam('verbose','yes')
		 
			 # Set inner photpars
			 for m in range(len(imgs)):
			   set_phot_params(fwhm[m],newdir)
			   if Ap1 == Ap2:
				 NumAps=[fwhm[m]*float(Ap1)]
				 skyrad = [fwhm[m]*float(Ap3)]
			   else:
				 NumAps=fwhm[m]*a_range
	   
			   if variable_aperture=="no":
				 iraf.photpars.setParam('apertur',a_range[j]) # List of aperture radii in scale units
				 print '\n Aperture: ' + str(a_range[j]) + '\tOuter radius: ' + str(a_range[k]) + '\tThickness: '+str(anwidth[l])+'\n'
			   else:
				 iraf.photpars.setParam('apertur',NumAps[j]) # List of aperture radii in scale units
				 iraf.fitskypars.setParam('annulus',skyrad[k]) # Inner radius of sky annulus in scale units
				 print '\n Aperture: ' + str(NumAps[j]) + '\tOuter radius: ' + str(skyrad[k]) + '\tThickness: '+str(anwidth[l])+'\n'
			   
			   print "Doing: ",j+1,"/",len(a_range),"\t\t\t",k+1,"/",len(skyrad),"\t\t\t",l+1,"/",len(anwidth)
			   print round((float(m)/len(imgs))*100.,1),"% done"
			   iraf.photpars.saveParList(filename="phot-loop.par")
			   iraf.phot.saveParList(filename="photpars.pars")


			   iraf.daophot.phot(

			   ####	OBJ parameters	     ####
			   str(directory+'/'+par.data_sub_dir[0]+'/'+imgs[m]),
			   coords=ds9,
			   output=mag_dir+"/default",
			   datapar="data.par",
			   centerp="center.par",
			   fitskyp="fitsky.par",
			   photpar="phot-loop.par",
			   interac="no",
			   verify="no",
			   update="no",
			   verbose="no")
			   #################################

			   
			   print '___________________________________________________________________\n'      
			   # Dumping data from the .mag.1 files into text files
			 

			 iraf.daophot.pdump(	# Getting the position and flux of stars

			 ####	OBJ parameters	     ####
			 mag_dir+"/*.mag.1",			#
			 fields="XCEN,YCEN,FLUX,MAG,MERR,MSKY",#
			 expr="yes",			#
			 headers="no",			#
			 Stdout=newdir+'/xyf.list')		#
			 #################################

			 os.system( 'mv *par* ' +params  )
            
  print "Photometry done\n"
			 
  return
