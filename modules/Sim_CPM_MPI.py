from mpi4py import MPI
import numpy as np
import mrcfile
import time
from pyrosetta import *
from pyrosetta.toolbox import cleanATOM
import cupy as cp
from cupy import cuda
import os
import sys
import argparse
import pandas as pd

# Setting up MPI and pi constant
pi = np.pi
comm = MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()

# Call for arguments
if rank == 0:
	parser = argparse.ArgumentParser(description='This is a Cryo-EM map simulator program developed by the Kellogg lab. Ahead are several argument flags presented to define a simulation of your choosing. If parameters are not set, default parameters will simulate apoferritin model included utilizing files from the template set.')
	parser.add_argument('-pdb_file',default='apoF__Q__emd_20521.pdb',type=str,help='Provide the path to the PDB file to be simulated')
	parser.add_argument('-mrc_file',default='apoF_unsharpened.mrc',type=str,help='Provide the path to the mrc file to determine output box shape and apix ratio')
	parser.add_argument('-output', default='OutputSim.mrc', type=str,help='Option to define name of output simulated file')
	parser.add_argument('--rescale',default=False,action='store_true',help='If used, will create an additional sim file output that rescales mrc values between 0 and 1')
	parser.add_argument('--procpergpu',action='store',type=int,default=1,help='State how many processes should be run per GPU used. Default 1.')
	parser.add_argument('--multinode', action='store_true', default=False, help='For use of MPI along with a multi-node GPU cluster. Option utilized for script usage with SDSC Comet')
	parser.add_argument('--force_apix',default=0,type=float,help='Force apix value for simulation, disregarding input mrc info')
	parser.add_argument('--force_box',default=0,type=int,help='Force box size to certain dimension, disregarding input mrc info')
	parser.add_argument('--HAmp',default=1.0,type=float,help='Define scattering factor amplitude for hydrogen atoms')
	parser.add_argument('--Hstd',default=3.15,type=float,help='Define scattering factor STDEV for hydrogen atoms')
	parser.add_argument('--OCAmp',default=1.35,type=float,help='Define scattering factor amplitude for backbone carbonyl oxygen atoms')
	parser.add_argument('--OCstd',default=4.15,type=float,help='Define scattering factor STDEV for backbone carbonyl oxygen atoms')
	parser.add_argument('--HeavyAmp',default=1.2,type=float,help='Define scattering factor amplitude for general heavy atoms')
	parser.add_argument('--Heavystd',default=3.2,type=float,help='Define scattering factor STDEV for general heavy atoms')
	parser.add_argument('--allatom',action='store_true',default = False,help = 'Reads scat. factor table provided and associates atoms not specifically defined with general value. In TYPE column, write "all-atom')
	parser.add_argument('--readetable',action='store',default='',type=str,help= 'Reads input scat. factor table for Gaussian calculations, please use CSV, use columns TYPE,SFAMP,SFSTD')
	args = parser.parse_args()
else:
	args = None


# Reads in electron-scattering factor table, and changes atom scattering factor arguments for use later during Gaussian addition, preliminary for now
def SFRead(CsvFile):
	Table = pd.read_csv(CsvFile).reset_index(drop = True)
	if(allatom == True):
		for i in range(0,len(Table)):
			if(Table.at[i,'Type'] == 'all-atom'):
				args.HAmp = args.HeavyAmp = args.OCAmp = Table.at[i,'SFAMP']
				args.Hstd = args.Heavystd = args.OCstd = Table.at[i,'SFSTD']
	for i in range(0,len(Table)):
		if(Table.at[i,'Type'] == 'H'):
			args.HAmp = Table.at[i,'SFAMP']
			args.Hstd = Table.at[i,'SFSTD']
		if(Table.at[i,'Type'] == 'OCbb'):
			args.OCAmp = Table.at[i,'SFAMP']
			args.OCstd = Table.at[i,'SFSTD']


#Can use this to scale values of MRC between 0 and 1
def MRCScale(npdata,filename):
    Min = np.amin(npdata)
    Max = np.amax(npdata)
    print(Max-Min)
    CurrData = np.array(npdata)
    CurrData = npdata
    CurrData = CurrData - Min
    ReviData = CurrData/(Max-Min)
    NewF = filename[:-4] + '-SCALED.mrc'
    with mrcfile.new(NewF,overwrite= True) as mrc:
        mrc.set_data(np.float32(ReviData))
        mrc._set_voxel_size(apix,apix,apix)
        mrc.close()

# Gathers information on every atom (coordinates and element-type) from pose
def CoordSet(pose):
	AtomArray = []
	print(str(pose.total_residue()))
	for presi in range(1,pose.total_residue()+1):
		for pres_atj in range(1,pose.residue(presi).natoms()):
			CurrArr = []
			at_type = pose.residue(presi).atom_type(pres_atj).atom_type_name()
			coords = pose.residue(presi).xyz(pres_atj)
			CurrArr.append(at_type)
			CurrArr.append(coords[0])
			CurrArr.append(coords[1])
			CurrArr.append(coords[2])
			AtomArray.append(CurrArr)
	return AtomArray

# Sets up ii,jj,kk
def OutmapSetUp(mrc,apix):
	mrcsz = mrc.data.shape
	outmap = np.zeros(mrcsz)
	outmap = np.zeros((BoxS,BoxS,BoxS))
	mrcsz = outmap.shape
	ij = np.fft.fftfreq(mrcsz[1],((1)))
	ijn = ij

	ii = np.zeros(mrcsz)
	jj = np.zeros(mrcsz)
	kk = np.zeros(mrcsz)

	for i in range(0,mrcsz[1]):
		ii[i,:,:] = ijn[i]
	for i in range(0,mrcsz[1]):
		jj[:,i,:] = ijn[i]
	for i in range(0,mrcsz[1]):
		kk[:,:,i] = ijn[i]

	ii = np.fft.fftshift(np.asarray(ii))
	jj = np.fft.fftshift(np.asarray(jj))
	kk = np.fft.fftshift(np.asarray(kk))

	return outmap, ii,jj,kk

# Actual gaussian function used to add atom gaussians to output map
def GaussForm(AtomicData):
	OutputArray = cp.zeros((BoxS,BoxS,BoxS))
	OutputArray = cp.array(OutputArray, dtype = np.complex64)
	scalefac = float(apix)
	for atom in AtomicData:
		#t1 = time.time()
		if(atom[0][0] == 'H'):
			scattering_params = cp.array([args.HAmp,args.Hstd])
		elif(atom[0] == 'OCbb'):
			scattering_params = cp.array([args.OCAmp,args.OCstd])
		else:
			scattering_params = cp.array([args.HeavyAmp,args.Heavystd])
		scattering_params = (scattering_params/scalefac)
		coords = atom[1:]
		center = cp.array([cp.float(coords[0]/apix),coords[1]/apix,cp.float(coords[2]/apix)])
		s = cp.float(1/scattering_params[1])
		ampl = cp.float((1/cp.sqrt(cp.power(2*pi,3)))*(1/cp.power(s,3)))
		coords = None
		OutputArray += ((cp.float(scattering_params[0]) * cp.fft.ifftshift(ampl* cp.exp(-cp.power(pi,2)*(cp.power(ii,2)+cp.power(jj,2)+cp.power(kk,2))/(2*cp.power(s,2)) - ( (2*pi)*1j*(ii*center[0]+jj*center[1]+kk*center[2]) )) )))
		center = None
		#t2 = time.time()
		#print('Atom Addition Time: ' + str(t2-t1))
		ToAdd, Rem, ampl, s, center, coords, scattering_params,t1,t2 = None,None,None,None,None,None,None,None,None
	OutputArray = cp.asnumpy(OutputArray)
	return OutputArray

# Broadcast arguments to every GPU
args = comm.bcast(args, root=0)
mrc = mrcfile.open(args.mrc_file)

# Set apix and boxsize. Uses mrc input, but can be overridden with optional arguments
if(args.force_apix == 0):
	apix = np.float32(mrc.voxel_size['x'])
else:
	apix = args.force_apix

if(args.force_box == 0):
	BoxS = len(mrc.data)
else:
	BoxS = args.force_box

# Splits up processes to separate GPU's, multinode refers to use in COMET cluster
if(args.multinode == False):
	Dev = int(rank/(args.procpergpu))
else:
	RankID = rank
	RankID = comm.gather(RankID,root = 0)
	if rank == 0:
		CudList = os.environ['CUDA_VISIBLE_DEVICES'].split(',')
		RankID.sort()
		Use = 0
		CurrCud = 0
		CudDict = {}
		for i in range(0,len(RankID)):
			CudDict[RankID[i]] = CudList[CurrCud]
			Use += 1
			if(Use == args.procpergpu):
				CurrCud += 1
				Use = 0
			if(CurrCud == 4):
				CurrCud = 0
	else:
		Use = None
		CudDict = None
		CurrCud = None
		CudList = None

	CudDict = comm.bcast(CudDict,root = 0)
	Dev = int(CudDict.get(rank))

print('Rank ID of this Process: ' + str(rank))
print('This rank will be associated with GPU: ' + str(Dev))

comm.Barrier()

# Rank 0 process does pyrosetta work and gathers atom coordinates/type
if rank == 0:
	init()
	cleanATOM(args.pdb_file)
	CleanFName = args.pdb_file[:-4] + '.clean.pdb'
	pose = pose_from_pdb(CleanFName)
	AtomData = CoordSet(pose)
	print('Total atom amount: ' + str(len(AtomData)))
	chunks = [[] for _ in range(size)]
	for i, chunk in enumerate(AtomData):
		chunks[i%size].append(chunk)
else:
	AtomData = None
	chunks = None

# Separates atom sets to each process equally
AtomData = comm.scatter(chunks, root = 0)
comm.Barrier()

# Actual Gaussian formation
with cp.cuda.Device(Dev):
	outmap,ii,jj,kk = OutmapSetUp(mrc,apix)
	ii = cp.asarray(ii,dtype=np.complex64)
	jj = cp.asarray(jj,dtype=np.complex64)
	kk = cp.asarray(kk,dtype=np.complex64)
	SumArray = GaussForm(AtomData)
	comm.Barrier()

#  Sums up values from every process' SumArray into the totals array
if rank == 0:
	totals = np.zeros_like(SumArray)
else:
	totals = None

comm.Reduce(SumArray,totals, op=MPI.SUM,root = 0)

# Inverse Fourier transform of totals array to get real-space data, then prepares output data
if rank == 0:
	outmap = np.real(np.transpose(np.fft.ifftn(totals)))
	if(args.rescale == True):
		MRCScale(outmap,args.output)
	with mrcfile.new(args.output,overwrite=True) as mrc:
		mrc.set_data(np.float32((outmap)))
		mrc._set_voxel_size(apix,apix,apix)
	mrc.close()
