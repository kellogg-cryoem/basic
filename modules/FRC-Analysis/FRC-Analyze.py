#!/usr/bin/env python
# coding: utf-8

import sys
import basicEM as basicEM
from PIL import Image
from scipy.ndimage.interpolation import shift as shft
import scipy
import pyem.star as s
import mrcfile
from scipy.constants import physical_constants
import numpy as np
import Simulate_EM_images as simim
from matplotlib import pyplot as plt
import fourier_ring_corr as fourier_ring_corr
from scipy.interpolate import RegularGridInterpolator
import os.path as path
from mpi4py import MPI
import math
import skimage
import itertools
from matplotlib import rcParams
import argparse
import pandas as pd


comm = MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()

if rank == 0:
	parser = argparse.ArgumentParser(fromfile_prefix_chars='@',description='Script that calculates and averages FRC from set of particles.')
	parser.add_argument('-star_file',default='run_data-cut81.star',type=str,help='Provide the path to the particle star file')
	parser.add_argument('-mrc_reference',default='postprocess_masked.mrc',type=str,help='Provide the path to the mrc file used as reference for FRC')
	parser.add_argument('-output', default='NyquistPlot.png', type=str,help='Option to define name of output image plot file')
	parser.add_argument('-particle_amount',default=-1,type=int,help='How many particles should be used from the star file? Takes first N particles from set. If default -1 is used, entire particle content of star file is utilized.')
	parser.add_argument('-particle_path',default='alignaver3/alignaver',type=str,help='Define path to particles')
	parser.add_argument('--missing_rlnimgname',default=False,action='store_true',help='Call this if the star file is missing the rlnimagename column. This will assume you have a particle stack with the particles in order according to the star file. Check out EMPIAR10026/EMD2788 to understand why I made this. Use --submrcs to direct to particle stack.')
	parser.add_argument('--submrcs', default = 'particles.mrcs',type = str,help = 'Use this to direct to a particle stack that contains only the picked particles to be inspected.')
	parser.add_argument('--reso_axis',default = False,action='store_true',help = 'If used, x-axis of plot is in units of resolution, if not default of nyquist fraction is used.')
	parser.add_argument('--defxaxis',default = False,action = 'store_true',help = 'X-axis ticks are set to default values that range from 25 to 1 angstrom.')
	parser.add_argument('--force_apix',default = 0,type = float,help = 'If particle stack does not include apix, use this to force a value')
	parser.add_argument('--xaxis_high',default = 55,type = float,help = 'Set upper limit for x axis value if using resolution units')
	parser.add_argument('--xaxis_low',default = 0.8,type = float,help = 'Set lower limit for x axis value if using resolution units')
	parser.add_argument('--xtick_high', default = 26, type = float,help = 'Set upper limit for x axis tick values if using resolution units')
	parser.add_argument('--xtick_low', default = 0.7, type = float,help = 'Set lower limit for x axis tick values if using resolution units')
	parser.add_argument('--report_data',default = '',type = str, help = 'Use this to create an output file that includes resolution/FRC datapoints. CSV formatted.')
	parser.add_argument('--qualcheck', default = False, action = 'store_true', help = 'In development. Use this to separate individual particles based on noise and when each ptcl crosses a specified threshold.')
	parser.add_argument('--qualthreshold',default = 0.1, type = float,help = 'Set up FRC threshold for qualcheck option.')
	parser.add_argument('--qualstartreso', default = 12.5, type = float, help = 'For quality check sets up upper limit of values to check depending on reso range')
	parser.add_argument('--qualendreso', default = 1.0, type = float, help = 'For quality check sets up lower limit of values to check depending on reso range')
	parser.add_argument('--qualrepo', default = 'ParticleQuality.csv', help = 'Reports results of qualcheck option if used to a csv. Use this to specify output csv name.')
	args = parser.parse_args()
else:
	args = None

args = comm.bcast(args, root=0)

def ReverseArray(Array): 
    return [element for element in reversed(Array)] 

def ValCross(Array,Value):
	Counter = 0
	FirstVal = Array[0]
	if(FirstVal > Value):
		AboveVal = True
	else:
		AboveVal = False
	for i in range(1,len(Array)):
		CurrCheck  = (Array[i] > Value)
		if(CurrCheck != AboveVal):
			Counter += 1
			AboveVal = CurrCheck
	return Counter

def LastTime(XValArr,YValArr,Threshold):
	RevY = ReverseArray(YValArr)
	RevX = ReverseArray(XValArr)
	for i in range(0,len(RevY)):
		if(RevY[i] > Threshold):
			ValFind = i
			break
	return RevX[i]

def FileFind(atomnum,starDF,subii=0):
	if(args.missing_rlnimgname == False):
		ptcl_ii_name = star.iloc[AtomID[ptcli]].rlnImageName
		ii_ndx = int(ptcl_ii_name[:ptcl_ii_name.find('@')])
		fn = ptcl_ii_name[ptcl_ii_name.find('/'):]
		fn = fn[fn.rfind('/'):]
		mrcfilename =  args.particle_path+fn
	else:
		mrcfilename = args.particle_path + args.submrcs
		ii_ndx = subii
	return mrcfilename, ii_ndx


if rank == 0:
	star = s.parse_star(args.star_file)
	if(args.particle_amount == -1):
		args.particle_amount = len(star)
	reftest = args.mrc_reference
	lam = simim.wavelength(200.0)
	C_s = 2.7 * 1e7 #C_s given in mm, convert to angstrom
	astigmatism = 0
	amplitude_contrast = 0.1

else:
	star = None
	reftest = None
	lam = None
	C_s = None
	astigmatism = None
	amplitude_contrast = None

args.particle_amount = comm.bcast(args.particle_amount,root = 0)
star = comm.bcast(star,root = 0)
reftest = comm.bcast(reftest,root = 0)
lam = comm.bcast(lam,root =0)
C_s = comm.bcast(C_s,root =0)
astigmatism = comm.bcast(astigmatism,root = 0)
amplitude_contrast = comm.bcast(amplitude_contrast,root =0)

f,fsc,ndxs = [],[],[]

# 26063 was original amount, took out 91 because 81 was acting up
if rank == 0:
	print('Amount of particles being investigated: ' + str(args.particle_amount))
	AtomID = range(0,args.particle_amount)
	chunks = [[] for _ in range(size)]
	dataspl = int(len(AtomID)/size)
	Counter = 0
	Remainder = args.particle_amount%size
	for i in range(0,size):
		for x in range(0,dataspl):
			chunks[i].append(AtomID[Counter])
			Counter += 1
		if(i < Remainder):
			chunks[i].append(AtomID[Counter])
			Counter += 1

else:
	AtomID = None
	chunks = None

AtomID = comm.scatter(chunks,root = 0)
from scipy.spatial.transform import Rotation as aR

if rank == 0:
	ptcli = 0
	if(args.missing_rlnimgname == False):
		fname, ii_ndx = FileFind(AtomID[0],star)
	else:
		fname, ii_ndx = FileFind(AtomID[0],star,subii = AtomID[0])
	if path.exists(fname) == False:
		print('Cannot find!')
		print(fname)
	if path.exists(fname):
		ptcl_imgs = mrcfile.open(fname)
		apix = np.float32(ptcl_imgs.voxel_size['x'])
		boxsize = np.float32(ptcl_imgs.data.shape[1]) 
		if(args.force_apix != 0):
			apix = args.force_apix
else:
	apix = None
	boxsize = None

apix = comm.bcast(apix,root = 0)
boxsize = comm.bcast(boxsize, root = 0)

for ptcli in range(0,len(AtomID)):
	print('This is the current particle being examined: ' + str(ptcli))
	if(args.missing_rlnimgname == False):
		fname, ii_ndx = FileFind(AtomID[ptcli],star)
	else:
		fname, ii_ndx = FileFind(AtomID[ptcli],star,subii = AtomID[ptcli])
	if path.exists(fname) == False:
		print('Cannot find!')
		print(fname)
		break
	if path.exists(fname):
		ptcl_imgs = mrcfile.open(fname)
		Rrot = aR.from_euler('z',star.iloc[AtomID[ptcli]].rlnAngleRot,degrees=True)
		Rtilt = aR.from_euler('y',star.iloc[AtomID[ptcli]].rlnAngleTilt,degrees=True)
		Rpsi = aR.from_euler('z',star.iloc[AtomID[ptcli]].rlnAnglePsi,degrees=True)

		defocusU = star.iloc[AtomID[ptcli]].rlnDefocusU

		R = Rrot * Rtilt * Rpsi
		R.as_dcm()

		#simulate_image takes defocus in Angstrom, usually given in nm
		prj_img,prj_img_ctf,R = simim.simulate_image_direct(reftest, R.as_dcm(), lam, defocusU, C_s)
		#x-y translation
		trans = [star.iloc[AtomID[ptcli]].rlnOriginX,star.iloc[AtomID[ptcli]].rlnOriginY]
	
		c = prj_img.shape 


		orig_x = np.arange(0,c[0]) 
		orig_y = np.arange(0,c[1])
		x = np.arange(0+trans[0],c[0]+trans[0],1)
		y = np.arange(0+trans[1],c[1]+trans[1],1)

		#required for shifting image by x and y
		my_interpolator = RegularGridInterpolator(points = [orig_x,orig_y], values=prj_img_ctf, bounds_error=False, fill_value=0)

		new_x,new_y = np.meshgrid(x,y)
		#print(trans)
		#this is where translation is applied? EK?
		new_im = my_interpolator((new_y,new_x))
		plt.imshow(new_im*-1,cmap='gray')
		if(len(ptcl_imgs.data.shape) > 2):
			ptcl_img = ptcl_imgs.data[int(ii_ndx)-1,:,:]
		else:
			ptcl_img = ptcl_imgs.data
		f_ii,fsc_ii = fourier_ring_corr.FSC(ptcl_img, new_im,disp=0)
		f.append(f_ii)
		fsc.append(fsc_ii)
		ndxs.append(ptcli)

print('If you are running this script and it is currently hanging here with no new output printed, there is likely something wrong with the star file as it struggles to find a particle from your mrcs')
comm.Barrier()

f = comm.gather(f,root = 0)
fsc = comm.gather(fsc, root = 0)
ndxs = comm.gather(ndxs,root = 0)

if rank == 0:
	f = list(itertools.chain.from_iterable(f))
	fsc = list(itertools.chain.from_iterable(fsc))
	ndxs = list(itertools.chain.from_iterable(ndxs))
	ax = plt.gca()

	fw = open('fsctest.dat','w')
	for i in range(0,len(fsc)):
		fw.write(str(fsc[i])+ '\n')
	fw.close()

	if(args.reso_axis == True):
		pixel_radii = np.arange(len(f[0]))+1
		x_axis_inverse =  pixel_radii/(boxsize*apix)
		x_axis = 1/(x_axis_inverse)
		ax.set_xlim(args.xaxis_high,args.xaxis_low)
		plt.xscale('log')
		xlocs, xlabs = plt.xticks()
		tickset = [50,37.5,25,12.5,9,6,3,1.5,1.2,1]
		ticks = [x for x in tickset if x <= args.xtick_high and x >= args.xtick_low]
		plt.xticks(ticks,ticks)
		plt.xlabel('Resolution (1/Å)')
		plt.ylabel('Fourier Ring Correlation')
	else:
		x_axis = f[0]
		plt.xlabel('Fraction Nyquist')
		plt.ylabel('Fourier Ring Correlation')

	meanfsc = np.zeros(fsc[0].shape)

	if(args.qualcheck == True):
		print('Starting quality testing...')
		QualDF = pd.DataFrame()
		StartReso = args.qualstartreso
		EndReso = args.qualendreso
		for i in range(0,len(x_axis)):
			if(x_axis[i] < StartReso):
				Startii = i-1
				break
		for i in range(0,len(x_axis)):
			if(x_axis[i] < EndReso):
				Endii = i
				break

		for i in range(0,len(fsc)):
			CurrFSC = fsc[i]
			QualDF.at[i,'Ptcl'] = i
			#print('This is the current Startii: ' + str(Startii))
			#print('This is the current Endii: ' + str(Endii))
			Crosses = ValCross(CurrFSC[Startii:Endii], args.qualthreshold)
			ResoFind = LastTime(x_axis[Startii:Endii],CurrFSC[Startii:Endii],args.qualthreshold)
			QualDF.at[i,'Crosses'] = Crosses
			QualDF.at[i,'ResoFind'] = ResoFind

		QualDF.to_csv(args.qualrepo)	
			
	cropped_fsc = np.zeros((args.particle_amount,len(fsc[0])))
	for ii in range(0,args.particle_amount):
	    print(ii)
	    meanfsc = meanfsc + fsc[ii]
	    cropped_fsc[ii,:] = fsc[ii] 

	meanfsc = np.dot(meanfsc,(1.0/args.particle_amount) )
	stddev = np.std(cropped_fsc,axis=0)

	if (args.report_data != ''):
		RepoDF = pd.DataFrame()
		for i in range(0,len(meanfsc)):
			RepoDF.at[i,'Resolution'] = x_axis[i]
			RepoDF.at[i,'MeanFRC'] = meanfsc[i]
		RepoDF.to_csv(args.report_data)
	#print(np.mean(stddev[1,]))

	#rcParams['font.family'] = 'sans-serif'
	#rcParams['font.sans-serif'] = ['Arial']

	#if(args.defxaxis == True):
	#	xlocs,xlabs = plt.xticks()
	#	print(xlocs)

	plt.rcParams.update({'font.size':18})
	plt.rcParams.update({'font.weight':'bold'})
	plt.ylim(0,0.7)
	plt.plot(x_axis,meanfsc,color='red',linewidth=4)
	ax.set_aspect('auto')

	import matplotlib.lines as mlines
	def newline(p1, p2):
	    ax = plt.gca()
	    xmin, xmax = ax.get_xbound()

	    if(p2[0] == p1[0]):
	        xmin = xmax = p1[0]
	        ymin, ymax = ax.get_ybound()
	    else:
	        ymax = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmax-p1[0])
	        ymin = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmin-p1[0])

	    l = mlines.Line2D([xmin,xmax], [ymin,ymax],linestyle='--',color='black')
	    ax.add_line(l)
	    return l

	newline([0, 0.1],[1, 0.1])
	plt.savefig(args.output)
	#plt.fill_between(f[0],meanfsc,meanfsc-0.22,meanfsc+0.22,alpha=0.