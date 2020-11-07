from mpi4py import MPI
import numpy as np
import mrcfile
import time
from pyrosetta import *
from pyrosetta.toolbox import cleanATOM
import cupy as cp
import os
import sys

pi = np.pi
apix=0.5332
mempool = cp.get_default_memory_pool()
pinned_mempool = cp.get_default_pinned_memory_pool()
comm = MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()
BoxS = 256

def ExpMRCScale(npdata,filename):
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


def OutmapSetUp(pose,mrc,apix):
	mrcsz = mrc.data.shape
	outmap = np.zeros(mrcsz)
	outmap = np.zeros((BoxS,BoxS,BoxS))
	mrcsz = outmap.shape
	ij = np.fft.fftfreq(mrcsz[1],((1/apix)))
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

    # scattering_params =  etbl.get(at_type)
	scattering_params = np.array([1,1]) #for now, use this

	return outmap, ii,jj,kk,scattering_params

def GaussForm(AtomicData):
	# At some point here, we need to add the portion that will find the scat factor based on atomtype
	# but for now...
	# scattering_params = etbl.get(at_type)
	#step = 0
	OutputArray = cp.zeros((BoxS,BoxS,BoxS))
	OutputArray = cp.array(OutputArray, dtype = np.complex64)
	scalefac = float(apix)
	for atom in AtomicData:
		#step += 1
		#t1 = time.time()
		if(atom[0][0] == 'H'):
			scattering_params = cp.array([1.0,3.15])
		elif(atom[0] == 'OCbb'):
			scattering_params = cp.array([1.35,4.15])
		else:
			scattering_params = cp.array([1.2,3.2])
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
		#print('Current atom: ' + str(step))
		ToAdd, Rem, ampl, s, center, coords, scattering_params,t1,t2 = None,None,None,None,None,None,None,None,None
	OutputArray = cp.asnumpy(OutputArray)
	#print('This is the size: ' + str(OutputArray.nbytes))
	return OutputArray

mrc = mrcfile.open('apoF-unsharpened.mrc')
init()
cleanATOM("apoF__Q__emd_20521.pdb")
pose = pose_from_pdb('apoF__Q__emd_20521.clean.pdb')


print('Rank ID of this Process: ' + str(rank))

Dev = int(rank+7)

if rank == 0:
	AtomData = CoordSet(pose)
	print('Total atom amount: ' + str(len(AtomData)))
	chunks = [[] for _ in range(size)]
	for i, chunk in enumerate(AtomData):
		chunks[i%size].append(chunk)
else:
	AtomData = None
	chunks = None

AtomData = comm.scatter(chunks, root = 0)
print('This rank will be associated with GPU: ' + str(Dev))
comm.Barrier()

with cp.cuda.Device(Dev):
	outmap,ii,jj,kk,scattering_params = OutmapSetUp(pose,mrc,apix)
	ii = cp.asarray(ii,dtype=np.complex64)
	jj = cp.asarray(jj,dtype=np.complex64)
	kk = cp.asarray(kk,dtype=np.complex64)
	SumArray = GaussForm(AtomData)
	comm.Barrier()

if rank == 0:
	totals = np.zeros_like(SumArray)
else:
	totals = None

comm.Reduce(SumArray,totals, op=MPI.SUM,root = 0)

if rank == 0:
	outmap = np.real(np.transpose(np.fft.ifftn(totals)))
	ExpMRCScale(outmap,'apoF-US-OGPDB-APIX05332-SF_FREQRevise.mrc')
	print(outmap)
	NewFName = 'apoF-US-OGPDB-APIX05332-SF_FREQRevise.mrc'
	with mrcfile.new(NewFName,overwrite=True) as mrc:
		mrc.set_data(np.float32((outmap)))
		mrc._set_voxel_size(apix,apix,apix)
	mrc.close()

