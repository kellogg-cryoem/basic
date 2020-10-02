from mpi4py import MPI
import numpy as np
import mrcfile
import time
from pyrosetta import *
from pyrosetta.toolbox import cleanATOM
import cupy as cp
import os
import sys

# This version of map simulator currently requires you to input values for hydrogen scattering factor parameters. As of 10/1/20, the optimized value is 0.75, 0.375.
# Declaring constants, receiving and formatting input
ScatParams = [float(i) for i in sys.argv[1:]]
pi = np.pi
apix=0.4

# Next two lines were used to determine memory usage, 3 lines after that were set up for further MPI-based commands
mempool = cp.get_default_memory_pool()
pinned_mempool = cp.get_default_pinned_memory_pool()
comm = MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()

# Used to convert PDB structure into set of atomic coordinates. It is turned into an array that can be scattered by MPI.
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

# Used for sampling distribution
def OutmapSetUp(pose,mrc,apix):
	mrcsz = mrc.data.shape
	outmap = np.zeros(mrcsz)
	
	ij = np.fft.fftfreq(mrcsz[1],1)
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

# Major function that takes up the most time, adds Gaussian per atom
def GaussForm(AtomicData,HParams):
	OutputArray = cp.zeros((340,340,340))
	OutputArray = cp.array(OutputArray, dtype = np.complex64)
	for atom in AtomicData:
		if(atom[0][0] == 'H'):
			scattering_params = cp.array(HParams)
		else:
			scattering_params = cp.array([1,3])
		coords = atom[1:]
		center = cp.array([cp.float(coords[0]/apix),coords[1]/apix,cp.float(coords[2]/apix)])
		s = cp.float(1/scattering_params[1])
		ampl = cp.float((1/cp.sqrt(cp.power(2*pi,3)))*(1/cp.power(s,3)))
		coords = None
		OutputArray += ((cp.float(scattering_params[0]) * cp.fft.ifftshift(ampl* cp.exp(-cp.power(pi,2)*(cp.power(ii,2)+cp.power(jj,2)+cp.power(kk,2))/(2*cp.power(s,2)) - ( (2*pi)*1j*(ii*center[0]+jj*center[1]+kk*center[2]) )) )))
		center = None
		ToAdd, Rem, ampl, s, center, coords, scattering_params,t1,t2 = None,None,None,None,None,None,None,None,None
	OutputArray = cp.asnumpy(OutputArray)
	return OutputArray

# As the first benchmark, used apoF model (1.2 A res) from Scheres et. al. See https://doi.org/10.1101/2020.05.22.110189 
mrc = mrcfile.open('apoF.mrc')
# Initializes pyrosetta
init()
cleanATOM("apoF__Q__emd_20521.pdb")
pose = pose_from_pdb('apoF__Q__emd_20521.clean.pdb')

print('Rank ID of this Process: ' + str(rank))

# Currently tested using RTX 2080Ti GPU's... each GPU is able to handle 2 processes of GaussForm at a time, please adjust the next line depending on your GPU specifications
Dev = int((rank/2)+1)

# Prepares coordiante data for scattering
if rank == 0:
	AtomData = CoordSet(pose)
	print('Current parameters being used for H: ' + str(ScatParams))
	print('Total atom amount: ' + str(len(AtomData)))
	chunks = [[] for _ in range(size)]
	for i, chunk in enumerate(AtomData):
		chunks[i%size].append(chunk)
else:
	AtomData = None
	chunks = None

AtomData = comm.scatter(chunks, root = 0)

# Tells you which process is associated with what GPU.
print('This rank will be associated with GPU: ' + str(Dev))

with cp.cuda.Device(Dev):
	outmap,ii,jj,kk= OutmapSetUp(pose,mrc,apix)
	ii = cp.asarray(ii,dtype=np.complex64)
	jj = cp.asarray(jj,dtype=np.complex64)
	kk = cp.asarray(kk,dtype=np.complex64)
	SumArray = GaussForm(AtomData,ScatParams)
	comm.Barrier()
	# Must keep barrier, otherwise code will crash if process finishes earlier compared to others

if rank == 0:
	totals = np.zeros_like(SumArray)
else:
	totals = None

# Returns and sums values from all GPU GaussForm functions
comm.Reduce(SumArray,totals, op=MPI.SUM,root = 0)

# Writes output
if rank == 0:
	outmap = np.real(np.transpose(np.fft.ifftn(totals)))
	print(outmap)
	NewFName = 'apoF-MPI-HMatrix-Heavy_1_3-HydA' + str(ScatParams[0]).replace('.','_').strip() + '-S' + str(ScatParams[1]).replace('.','_').strip() + '.mrc'
	with mrcfile.new(NewFName,overwrite=True) as mrc:
		mrc.set_data(np.float16((outmap)))
		# Sets header information for output mrc file
		mrc._set_voxel_size(apix,apix,apix)
	mrc.close()
