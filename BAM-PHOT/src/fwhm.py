import os.path
import numpy as np
import param as par
import sys
from pyraf import iraf

def fwhm(directory):
  iraf.noao(_doprint=0)
  iraf.obsutil(_doprint=0)
  iraf.unlearn("psfmeasure")
  
  rootdir = par.rootdir[0]
  
  dir_contents = os.listdir(directory+'/'+par.data_sub_dir[0])
  ds9s = [fn for fn in dir_contents if fn.startswith(par.ds9_name[0]) and fn.endswith(par.ds9_name[1])]
  ds9 = directory+'/'+par.data_sub_dir[0]+'/'+ds9s[0]

  print "\nMeasuring the fwhm of selected reference stars..."
  iraf.noao.obsutil.psfmeasure.coords = "mark1"
  iraf.noao.obsutil.psfmeasure.display = "no"
  iraf.noao.obsutil.psfmeasure.size = "FWHM"
  iraf.noao.obsutil.psfmeasure.scale = "1"
  iraf.noao.obsutil.psfmeasure.radius = "11"
  iraf.noao.obsutil.psfmeasure.swidth = "15"
  iraf.noao.obsutil.psfmeasure.iterati = "1"
  iraf.noao.obsutil.psfmeasure.imagecu = ds9
  iraf.noao.obsutil.psfmeasure.graphcur = rootdir+"/src/q.reg"
  iraf.noao.obsutil.psfmeasure.logfile = directory+"/"+par.data_sub_dir[0]+"/fwhm.log"

  if os.path.isfile(directory+"/"+par.data_sub_dir[0]+"/fwhm.log"):
    os.system("rm "+directory+"/"+par.data_sub_dir[0]+"/fwhm.log")
  
  imgs = [fn for fn in dir_contents if fn.startswith(par.image_name[0]) and fn.endswith(par.image_name[1])]
  imgs = sorted(imgs)

  for i in range(len(imgs)):
    iraf.noao.obsutil.psfmeasure.imagecu = ds9
    iraf.noao.obsutil.psfmeasure(directory+'/'+par.data_sub_dir[0]+'/'+imgs[i])
  
  N_stars = len(ds9)
  fwhm = [[] for _ in range(len(imgs))]
  values = [ line for line in open(directory+"/"+par.data_sub_dir[0]+"/fwhm.log") if '.' in line and 'NOAO' not in line]
  j = 0
  for i in range(len(values)):
    if values[i][2:9] == 'Average':
      j += 1
    if values[i][2:9] != 'Average': 
      fwhm[j].append(float(values[i][41:47]))

  f = open(directory+'/'+par.data_sub_dir[0]+'/fwhm.txt', 'w+')
  for i in range(len(imgs)):
    fwhm[i] = np.array(fwhm[i])
    print >> f, np.median(fwhm[i]), imgs[i]
  f.close()
