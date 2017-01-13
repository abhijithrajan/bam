#! /usr/bin/env python
pixel_scale = [0.288]
read_noise = [11.34]
gain = [5.385]
JD_header = ["MJD-OBS"]
airmass = ["HIERARCH ESO TEL AIRM START"]
rootdir = ["/Volumes/BigBackup/MMT/NTT_phot"] # i.e "/home/user/data" Basically the path where the paw.py file is.
datadir = ["data"]
binned_data = [False]
data_sub_dir = ["data_unbinned"]
ds9_name = 		["ds9","v2.reg"] # First and last parts of the coordinate 'ds9.reg' file.
image_name = 		["2M1",".fits"]	# First and last parts of the .fits files.
variable_aperture = ["yes"]
airmass_detrend = 	["yes"]

std_cutoff = [0.0]
cut_offs = [2.0,0.0,10,2,50.,10000.]	# [?,MAD cut_off (Not used),number of starting refs,min number of refs,min pix value,max pix value]

ignore_objects = []#34,38,51,56,60,63,71,72,77]
