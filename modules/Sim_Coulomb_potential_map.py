#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cupy
import mrcfile
import numpy


# In[2]:


def read_etable(etbl_file):
    f = open('examples/etable.txt')
    ln = f.readlines()
    etbl = {}
    for ii in range(0,len(ln)):
        itms = ln[ii].split()
        if( (len(itms) >= 3) & (itms[0] != "#") & (itms[0] != "NAME") & (itms[0] != "###")):
            #print(itms[0]+" "+itms[1]+" "+itms[2]+"\n")
            etbl[itms[0]] = [float(itms[2]),float(itms[3])]
    return etbl

ee = read_etable('examples/etable.txt')
arr=ee.get('CNH2')
print(arr)


# In[3]:


def sim_map_lite(pose,ii,jj,kk,boxsize,apix,etbl):
    #fourier indices as input(ii,jj,kk)
    import cupy
    outmap = cupy.zeros((boxsize,boxsize,boxsize))
    
    ii = cupy.asarray(ii)
    jj = cupy.asarray(jj)
    kk = cupy.asarray(kk)
    
    pi = cupy.pi

    for presi in range(1,pose.total_residue()):
    #get number of atoms
    #print(presi)
        for pres_atj in range(1,pose.residue(presi).natoms()):
            #for each atom, simulate gaussian with appropriate electron scattering factor
            #at_type = pose.residue(presi).atom_type(pres_atj)
            #print('this is the atom-type:' + at_type.str())
        
            #scattering_params = etbl.get(at_type)
            scattering_params = cupy.array([1, 1])
            coords = pose.residue(presi).xyz(pres_atj) 
            center = cupy.array( [ coords[0]/apix, coords[1]/apix, coords[2]/apix ])
            #print(center)
            s = cupy.float(1/scattering_params[1])

            ampl = cupy.float((1/cupy.sqrt(cupy.power(2*pi,3)))*(1/cupy.power(s,3)))
        
            ii_ampl = cupy.float(scattering_params[0])
            
            outmap = outmap +  ii_ampl * cupy.fft.ifftn(cupy.fft.ifftshift(ampl* cupy.exp(-cupy.power(pi,2)*(cupy.power(ii,2)+cupy.power(jj,2)+cupy.power(kk,2))/(2*cupy.power(s,2)) - ( (2*pi)*1j*(ii*center[0]+jj*center[1]+kk*center[2]) )) ))
            
    outmap = numpy.real(cupy.asnumpy(cupy.transpose(outmap)))
    return outmap


# In[4]:


def sim_map(pose,mrc,apix,etbl):

    #import mrcfile as mrc
    #from pyrosetta import *
    #init()


    #import numpy
    import cupy
    #from pyrosetta.toolbox import cleanATOM

    #mrc = mrcfile.open('examples/1ye3.mrc')
    #cleanATOM('examples/1ye3.pdb')
    #pose = pose_from_pdb('examples/1ye3.clean.pdb')

    #etbl = read_etable('examples/etable.txt')
    #apix=1

    mrcsz = mrc.data.shape
    print(mrcsz)
    outmap = cupy.zeros(mrcsz)
    
    ij = cupy.fft.fftfreq(mrcsz[1],1)
    
    ii = numpy.zeros(mrcsz)
    jj = numpy.zeros(mrcsz)
    kk = numpy.zeros(mrcsz)
    
    for ii_ndx in range(0,mrcsz[1]):
        for jj_ndx in range(0,mrcsz[1]):
            for kk_ndx in range(0,mrcsz[1]):
                ii[ii_ndx,jj_ndx,kk_ndx] = ij[ii_ndx]
                jj[ii_ndx,jj_ndx,kk_ndx] = ij[jj_ndx]
                kk[ii_ndx,jj_ndx,kk_ndx] = ij[kk_ndx]

    print('here')
    ii = cupy.fft.fftshift(cupy.asarray(ii))
    jj = cupy.fft.fftshift(cupy.asarray(jj))
    kk = cupy.fft.fftshift(cupy.asarray(kk))

    pi = cupy.pi

    for presi in range(1,pose.total_residue()):
    #get number of atoms
    #print(presi)
        for pres_atj in range(1,pose.residue(presi).natoms()):
            #for each atom, simulate gaussian with appropriate electron scattering factor
            #at_type = pose.residue(presi).atom_type(pres_atj)
            #print('this is the atom-type:' + at_type.str())
        
            #scattering_params = etbl.get(at_type)
            scattering_params = cupy.array([1, 1])
            coords = pose.residue(presi).xyz(pres_atj) 
            center = cupy.array( [ coords[0]/apix, coords[1]/apix, coords[2]/apix ])
            #print(center)
            s = cupy.float(1/scattering_params[1])

            ampl = cupy.float((1/cupy.sqrt(cupy.power(2*pi,3)))*(1/cupy.power(s,3)))
        
            ii_ampl = cupy.float(scattering_params[0])
            
            outmap = outmap +  ii_ampl * cupy.fft.ifftn(cupy.fft.ifftshift(ampl* cupy.exp(-cupy.power(pi,2)*(cupy.power(ii,2)+cupy.power(jj,2)+cupy.power(kk,2))/(2*cupy.power(s,2)) - ( (2*pi)*1j*(ii*center[0]+jj*center[1]+kk*center[2]) )) ))
            
    outmap = numpy.real(cupy.asnumpy(cupy.transpose(outmap)))
    return outmap


# In[5]:


from matplotlib import pyplot as plt
inmrc = mrcfile.open('examples/1ye3.mrc')
plt.imshow(inmrc.data[64,:,:])

mrc = mrcfile.open('examples/1ye3.mrc')
from pyrosetta import *
import numpy as np
init()
from pyrosetta.toolbox import cleanATOM
cleanATOM("./examples/1ye3.pdb")
pose = pose_from_pdb("./examples/1ye3.clean.pdb")

etbl = read_etable('examples/etable.txt')
apix=1

omap = sim_map(pose,mrc,apix,etbl)


# In[6]:


from matplotlib import pyplot as plt

with mrcfile.new('outputtest.mrc',overwrite=True) as mrc:
    mrc.set_data(numpy.float16((omap)))
mrc.close()

plt.imshow(((omap[64,:,:])))






