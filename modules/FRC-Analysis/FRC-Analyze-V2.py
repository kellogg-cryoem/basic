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
import fourier_ring_corr
from mpi4py import MPI
import math
import skimage
import itertools
from matplotlib import rcParams
import argparse
import pandas as pd
import CTF_calculator as ctfcalc
from scipy.io import savemat
import cv2


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
	parser.add_argument('--projected',default = False,action='store_true', help = 'Set defocus to 0, dummy value for occassions like when STAR file does not include it, such as when using relion_project images.')
	parser.add_argument('--starapix', default = 0, type = float, help = 'If parse_star is giving you an imagepixel error, use this to force it to a value.')
	parser.add_argument('--plotprefix',default = '',type = str,help = 'If option is specified, plots for CTF, RefParticles, and SimParticle. Use mpirun -n 1 to produce all particle images.')
	parser.add_argument('--refptclresize', default = False, action = 'store_true',help = 'During ctf application, ctf is expanded to 3x of original size, with reference image subjected to same method before convolution.')
	parser.add_argument('--ptclimgpad', default = False, action = 'store_true', help = 'Variation of CTF application. Pads ptcl img in real space to counter aliasing effects of CTF with small boxsize')
	parser.add_argument('--ptclimgpad2', default = False, action = 'store_true', help = 'Variation of CTF application. Pads ptcl img in real space to counter aliasing effects of CTF with small boxsize')
	parser.add_argument('--apixctf', default = False, action = 'store_true', help = 'If used, the apix value inputted to simulate CTF is the same as the mrc apix, and is not changed to remove artifacts')
	parser.add_argument('--particle_mrcs',default = '', type = str, help = 'Adds reference particle images together in mrcs.')
	parser.add_argument('--outputtxt', default = '', type = str, help = 'Save FRC output as text file')
	parser.add_argument('--sim_mrcs',default = '', type = str, help = ' Adds simulated images into mrcs file')
	parser.add_argument('--savestar', default = '', type = str, help = 'Save star file, converted to DF, as a csv.')
	parser.add_argument('--simload', default= '',type = str, help = 'If you got simulated images alreacy created, use this to load the mrcs.')
	parser.add_argument('--mask', default = '', type = str, help = 'Use option to add mask to raw image. For now, can only be used with simload option.')
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
	if(args.savestar != ''):
		star.to_csv(args.savestar)
	if(args.particle_amount == -1):
		args.particle_amount = len(star)
	reftest = args.mrc_reference
	lam = simim.wavelength(300.0)
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

if rank == 0:
	print('Amount of particles being investigated: ' + str(args.particle_amount))
	AtomID = range(0,args.particle_amount)
	chunks = [[] for _ in range(size)]
	dataspl = int(len(AtomID)/size)
	Counter = 0
	Remainder = args.particle_amount%size

	if(args.particle_mrcs != ''):
		NewMRC = mrcfile.new(args.particle_mrcs, overwrite = True)
		WriteMRC = []
	if(args.sim_mrcs != ''):
		SimMRCS = mrcfile.new(args.sim_mrcs,overwrite = True)
		WriteSim = []

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
	NewMRC = None
	WriteMRC = None

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
	if(args.mask != ''):
		mask = mrcfile.open(args.mask)
		mask = mask.data
	if path.exists(fname):
		ptcl_imgs = mrcfile.open(fname)
		apix = np.float32(ptcl_imgs.voxel_size['x'])
		boxsize = np.float32(ptcl_imgs.data.shape[1]) 
		if(args.force_apix != 0):
			apix = args.force_apix
else:
	apix = None
	boxsize = None
	if(args.mask != ''):
		mask = None

apix = comm.bcast(apix,root = 0)
boxsize = comm.bcast(boxsize, root = 0)
if(args.mask != ''):
	mask = comm.bcast(mask, root = 0)

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

		if(args.projected == False):		
			defocusU = star.iloc[AtomID[ptcli]].rlnDefocusU
			defocusV = star.iloc[AtomID[ptcli]].rlnDefocusV
			defocusang = star.iloc[AtomID[ptcli]].rlnDefocusAngle
		else:
			defocusU = 0
			defocusV = 0
			defocusang = 0

		R = Rrot * Rtilt * Rpsi
		R.as_dcm()


		if(args.simload != ''):

			if(len(ptcl_imgs.data.shape) > 2):
				try:
					ptcl_img = ptcl_imgs.data[int(ii_ndx)-1,:,:]
				except:
					print('Ptcl_img not found.')
					continue
			else:
				ptcl_img = ptcl_imgs.data

			if(args.mask != ''):
					#ptcl_fft = np.fft.fftn(ptcl_img)
					#ptcl_fft_shift = np.fft.fftshift(ptcl_fft)
					ptcl_mask,ptcl_mask_ctf,R,OGctf = simim.simulate_image_direct(args.mask, R.as_dcm(), lam,defocusU,defocusV,defocusang, C_s)
					#ptcl_mask_fft = np.fft.fftn(ptcl_mask)
					#ptcl_mask_fft_shift = np.fft.fftshift(ptcl_mask_fft)
					#ptcl_fft_ptcl_mask_multi = ptcl_fft_shift * ptcl_mask_fft_shift
					ptcl_multiply = ptcl_img * ptcl_mask
					#ptcl_masked = np.fft.ifftn(ptcl_multiply_ifftshift)
					ptcl_img = ptcl_multiply

			sim_imgs = mrcfile.open(args.simload)
			if(len(sim_imgs.data.shape) > 2):
				try:

					ifft_ctf = sim_imgs.data[AtomID[ptcli],:,:]
				except:
					print('Image not found.')
					print('This is the ii_ndx looked for: ' + str(ii_ndx))
					continue
			else:
				ifft_ctf = sim_imgs.data

			if (args.plotprefix != ''):
				PtclFile = args.plotprefix + '-RefPtcl-' + str(AtomID[ptcli]) + '.png'
				SimFile = args.plotprefix + '-SimPtcl-' + str(AtomID[ptcli]) + '.png'

				plt.imsave(PtclFile,(np.real(ptcl_img).astype(np.float32)))
				plt.imsave(SimFile,ifft_ctf)


			f_ii,fsc_ii = fourier_ring_corr.FSC(ptcl_img, ifft_ctf,disp=0)
			f.append(f_ii)
			fsc.append(fsc_ii)
			ndxs.append(ptcli)
			continue


		#simulate_image takes defocus in Angstrom, usually given in nm
		prj_img,prj_img_ctf,R,OGctf = simim.simulate_image_direct(reftest, R.as_dcm(), lam,defocusU,defocusV,defocusang, C_s)
		#print('defocusU: ' + str(defocusU))
		#print('defocusV: ' + str(defocusV))
		#print('defocusang: ' + str(defocusang))
		#print('Cs: ' + str(C_s))
		#x-y translation
		trans = [star.iloc[AtomID[ptcli]].rlnOriginX,star.iloc[AtomID[ptcli]].rlnOriginY]
		c = prj_img.shape 
		

		orig_x = np.arange(0,c[0]) 
		orig_y = np.arange(0,c[1])
		x = np.arange(0+trans[0],c[0]+trans[0],1)
		y = np.arange(0+trans[1],c[1]+trans[1],1)

		#required for shifting image by x and y
		#my_interpolator = RegularGridInterpolator(points = [orig_x,orig_y], values=prj_img_ctf, bounds_error=False, fill_value=0)
		#my_interpolator = RegularGridInterpolator(points = [orig_x,orig_y], values=prj_img_ctf, bounds_error=False, fill_value=0)
		my_interpolator = RegularGridInterpolator(points = [orig_x,orig_y], values=prj_img, bounds_error=False, fill_value=0)
		my_interpolator = RegularGridInterpolator(points = [orig_x,orig_y], values=prj_img, bounds_error=False, fill_value=0)
		

		new_x,new_y = np.meshgrid(x,y)

		#this is where translation is applied? EK?
		new_im = my_interpolator((new_y,new_x))
		#plt.imshow(new_im*-1,cmap='gray')
		if(len(ptcl_imgs.data.shape) > 2):
			try:
				ptcl_img = ptcl_imgs.data[int(ii_ndx)-1,:,:]
			except:
				continue
		else:
			ptcl_img = ptcl_imgs.data

		angpix = np.float32(ptcl_imgs.voxel_size['x'])
		if angpix == 0:
			angpix = args.force_apix
		ogdim = new_im.shape[0]
		ctf_samples = ctfcalc.get_ctf(angpix,ogdim,defocusU,defocusV,defocusang,C_s,300,0.1,0)
		prj_fft = np.fft.fftn(new_im)
		#prj_fft = np.pad(prj_fft,ogdim/2)
		prj_fft = np.fft.fftshift(prj_fft)
		#prj_fft = np.pad(prj_fft,int(ogdim))
		prj_fft_ctf = prj_fft *ctf_samples
		#prj_fft_ctf = prj_fft_ctf[int(ogdim):int(-ogdim),int(ogdim):int(-ogdim)]
		prj_fft_ctf = np.fft.ifftshift(prj_fft_ctf)
		ifft_ctf = np.real(np.fft.ifftn((prj_fft_ctf)))
		#ifft_ctf = ifft_ctf[ogdim/2:-ogdim/2,ogdim/2:-ogdim/2]

		if(args.refptclresize == True):
			newdim = ogdim * 3
			ctf_samples = ctfcalc.get_ctf(angpix,newdim,defocusU,defocusV,defocusang,C_s,300,0.1,0)
			new_im = cv2.resize(new_im,dsize = (newdim,newdim))
			prj_fft = np.fft.fftn(new_im)
			prj_fft = np.fft.fftshift(prj_fft)
			prj_fft_ctf = prj_fft * ctf_samples
			prj_fft_ctf = np.fft.ifftshift(prj_fft_ctf)
			ifft_ctf = np.real(np.fft.ifftn((prj_fft_ctf)))
			ifft_ctf = cv2.resize(ifft_ctf,dsize = (ogdim,ogdim))

			CurrCTF = args.plotprefix + '-CTF-' + str(ptcli) + '.png'
			PtclFile = args.plotprefix + '-RefPtcl-' + str(ptcli) + '.png'
			SimFile = args.plotprefix + '-SimPtcl-' + str(ptcli) + '.png'

			plt.imsave(CurrCTF,ctf_samples)
			plt.imsave(PtclFile,ptcl_img)
			plt.imsave(SimFile,ifft_ctf)

		if(args.ptclimgpad == True):
			newdim = ogdim*3
			ctf_samples = ctfcalc.get_ctf(angpix,newdim,defocusU,defocusV,defocusang,C_s,300,0.1,0)
			new_im = np.pad(new_im,ogdim)
			print(new_im.shape)
			print(ctf_samples.shape)
			prj_fft = np.fft.fftn(new_im)
			prj_fft = np.fft.fftshift(prj_fft)
			prj_fft_ctf = prj_fft * ctf_samples
			prj_fft_ctf = np.fft.ifftshift(prj_fft_ctf)
			ifft_ctf = np.real(np.fft.ifftn((prj_fft_ctf)))
			#ifft_ctf = ifft_ctf[int(ogdim):int(-ogdim),int(ogdim):int(-ogdim)]

			CurrCTF = args.plotprefix + '-CTF-' + str(ptcli) + '.png'
			PtclFile = args.plotprefix + '-RefPtcl-' + str(ptcli) + '.png'
			SimFile = args.plotprefix + '-SimPtcl-' + str(ptcli) + '.png'

			plt.imsave(CurrCTF,ctf_samples)
			plt.imsave(PtclFile,ptcl_img)
			plt.imsave(SimFile,ifft_ctf)

		if(args.ptclimgpad2 == True):
			newdim = ogdim*3
			ctf_samples = ctfcalc.get_ctf(angpix,newdim,defocusU,defocusV,defocusang,C_s,300,0.1,0)
			#new_im = np.pad(new_im,ogdim)
			save_im = new_im
			prj_fft = np.fft.fftn(new_im)
			prj_fft = np.fft.fftshift(prj_fft)
			prj_fft = np.pad(prj_fft,ogdim)
			prj_fft_ctf = prj_fft * ctf_samples
			prj_fft_ctf = np.fft.ifftshift(prj_fft_ctf)
			ifft_ctf = np.real(np.fft.ifftn((prj_fft_ctf)))
			ifft_ctf = ifft_ctf[int(ogdim):int(-ogdim),int(ogdim):int(-ogdim)]

			CurrCTF = args.plotprefix + '-CTF-' + str(ptcli) + '.png'
			PtclFile = args.plotprefix + '-RefPtcl-' + str(ptcli) + '.png'
			SimFile = args.plotprefix + '-SimPtcl-' + str(ptcli) + '.png'

			plt.imsave(CurrCTF,ctf_samples)
			plt.imsave(PtclFile,ptcl_img)
			plt.imsave(SimFile,ifft_ctf)

		if(args.particle_mrcs != ''):
			if(rank == 0):
				WriteMRC.append(ptcl_img)
		if(args.sim_mrcs != ''):
			if(rank == 0):
				WriteSim.append(ifft_ctf)

		f_ii,fsc_ii = fourier_ring_corr.FSC(ptcl_img, ifft_ctf,disp=0)
		ogdim = new_im.shape[0]

		#new_im = np.pad(new_im,new_im.shape[0])
		#new_im = new_im[ogdim:-ogdim,ogdim:-ogdim]
		if (args.plotprefix != '' and args.ptclimgpad2 != True and args.refptclresize != True and args.ptclimgpad != True):
			CurrCTF = args.plotprefix + '-CTF-' + str(ptcli) + '.png'
			PtclFile = args.plotprefix + '-RefPtcl-' + str(ptcli) + '.png'
			SimFile = args.plotprefix + '-SimPtcl-' + str(ptcli) + '.png'

			CTFArr = args.plotprefix + '-CTFArr-' + str(AtomID[ptcli]) + '.mat'
			SimArr = args.plotprefix + '-SimArr-' + str(AtomID[ptcli]) + '.mat'
			RawArr = args.plotprefix + '-RawArr-' + str(AtomID[ptcli]) + '.mat'
			RawSimArr = args.plotprefix + '-RawSimArr-' + str(AtomID[ptcli]) + '.mat'

			savemat(RawSimArr,mdict={'arr':new_im})

			# savemat(CTFArr,mdict={'arr':ctf_samples})
			# savemat(SimArr,mdict={'arr':ifft_ctf})

			# plt.imsave(CurrCTF,ctf_samples)
			# plt.imsave(PtclFile,ptcl_img)
			# plt.imsave(SimFile,ifft_ctf)


		# Used 4 next commands to print 
		#CTFMat = 'Ptcl9-Test1-CTF-' + str(ptcli) + '-Array.mat'
		#SimMat = 'Ptcl9-Test1-SimCTFPtcl-' + str(ptcli) + '-Array.mat'
		#savemat(CTFMat,mdict={'arr':ctf_samples})
		#savemat(SimMat,mdict={'arr':ifft_ctf})

		f.append(f_ii)
		fsc.append(fsc_ii)
		ndxs.append(ptcli)

print('If you are running this script and it is currently hanging here with no new output printed, there is likely something wrong with the star file as it struggles to find a particle from your mrcs')
comm.Barrier()

f = comm.gather(f,root = 0)
fsc = comm.gather(fsc, root = 0)
ndxs = comm.gather(ndxs,root = 0)

if rank == 0:

	if(args.particle_mrcs != ''):
		WriteMRC = np.array(WriteMRC)
		NewMRC.set_data(WriteMRC)
		NewMRC.close()

	if(args.sim_mrcs != ''):
		WriteSim = np.array(WriteSim)
		SimMRCS.set_data(WriteSim.astype(np.float32))
		SimMRCS.close()

	f = list(itertools.chain.from_iterable(f))
	fsc = list(itertools.chain.from_iterable(fsc))
	ndxs = list(itertools.chain.from_iterable(ndxs))
	ax = plt.gca()

	if(args.reso_axis == True):
		pixel_radii = np.arange(len(f[0]))+1
		x_axis_inverse =  pixel_radii/(boxsize*apix)
		x_axis = 1/(x_axis_inverse)
		ax.set_xlim(args.xaxis_high,args.xaxis_low)
		plt.xscale('log')
		xlocs, xlabs = plt.xticks()
		tickset = [50,37.5,25,12.5,9,6,3,1.5,1.2,1,0.8,0.6,0.5]
		ticks = [x for x in tickset if x <= args.xtick_high and x >= args.xtick_low]
		plt.xticks(ticks,ticks)
		plt.xlabel('Resolution (1/Å)')
		plt.ylabel('Fourier Ring Correlation')

		if(args.outputtxt != ''):
			XAx = args.outputtxt + '-Ptcl' + str() + '-XFreq.dat'
			f = open(XAx,'w')
			for val in x_axis:
				f.write(str(val) + '\n')
			f.close()
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
	for ii in range(0,len(fsc)):
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

	if(args.outputtxt != ''):
		#XAx = args.outputtxt + '-Ptcl' + str() + '-XFreq.dat'
		FSCfile = args.outputtxt + '-FSC.dat'
		#f = open(XAx,'w')
		#for val in x_axis:
		#	f.write(val + '\n')
		#f.close()
		np.savetxt(FSCfile,meanfsc)

	plt.rcParams.update({'font.size':18})
	plt.rcParams.update({'font.weight':'bold'})
	plt.ylim(0,1.1)
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