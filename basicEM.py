#!/usr/bin/env python
# coding: utf-8



import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import mrcfile as mrc




def square(m):
    #m is a numpy array
    import numpy as np
    import mrcfile 
    mindim = np.min(m.shape)
    return m[0:mindim,0:mindim]




def normalize(m):
    import numpy as np
    return np.dot(m,1/(sum(sum(m))))




def fft_crop(m,scale):
    #assumes im is a mrc object
    import numpy as np
    
    newsize = np.ndarray.astype(np.dot(m.data.shape,scale),int)
    #newsize = (np.dot(m.data.shape,0.1))
    c = np.ndarray.astype(newsize / 2,int)
    s = np.ndarray.astype((c - (newsize[0]/2)) +1,int)
    e = np.ndarray.astype((c + (newsize[1]/2)),int)
    fftim = np.fft.fftshift(np.fft.fft2(m.data))
    newim = np.real(np.fft.ifft2(np.fft.ifftshift(fftim[s[0]:e[0],s[1]:e[1]])))
    return newim



def fftim(m):
    import numpy as np
    return np.fft.fftshift(np.fft.fft2(m))
    



def gaussC( x,y,sigma,center):
    import numpy as np
    xc = center[0]
    yc = center[1]
    exponent = (np.divide((np.square(x-xc) + np.square(y-yc)),(2*sigma)))
    return np.exp(-exponent)

def gauss2d(sz,sigma,center):
    import numpy as np
    ndx = np.linspace(0,sz,sz)
    [R,C] = np.meshgrid(ndx,ndx);
    return gaussC(R,C,sigma,center)



def circmask(sz,radius,center):
    import numpy as np
    ndx = np.linspace(0,sz,sz)
    [X,Y] = np.meshgrid(ndx,ndx)
    mask = (np.square(X-center[0]) + np.square(Y-center[1])) <= radius**2
    mask = 1*mask.astype(float)
    return mask




def show_power_spectra(im,apix):
    get_ipython().run_line_magic('matplotlib', 'inline')
    from matplotlib import pyplot as plt
    import numpy as np
    plt.rcParams['figure.figsize'] = [25, 25]

    #1. FFT
    fft_im = ((fftim(square(im))))
    #2. crop FFT to 3 Angstrom
    ndx = np.dot( fft_im.shape, apix/3 )
    c = np.dot( fft_im.shape , 0.5 )
    s = c - ndx
    e = c + ndx
    s = np.ndarray.astype(s,int)
    e = np.ndarray.astype(e,int)
    fft_im = fft_im[ (s[0]):(e[0]), (s[1]):(e[1])]
    #3. take absolute value of FFT
    fft_im = np.abs(fft_im)
    #steps 4-6 are only for displaying the FFT
    #4. fft the FFT
    fft_squared = fftim(fft_im)
    #5. crop to 512
    c = np.ndarray.astype( np.dot( fft_squared.shape, 0.5 ), int )
    fft_squared = fft_squared[(c[0]-256):(c[0]+256),(c[1]-256):(c[1]+256)]
    #6. inverse the fft(fft) to get the fft
    TheOriginalfft = np.fft.ifft2(np.fft.ifftshift(fft_squared))
    sz = TheOriginalfft.shape
    center = np.dot(sz,0.5)
    #6b. high-pass filter the fft
    sigma = 100000
    hpmask = gauss2d(sz[0],sigma,center)
    hpmaskfft = np.multiply(abs(TheOriginalfft),hpmask)
    #7. mask out central peak
    #gmask=np.abs(np.add( np.dot( gauss2d(sz[0],10,center), -1 ), 1)) 
    cm=np.add( np.dot( circmask(sz[0],10,center), -1), 1)
    print(cm[256,256])
    TheOriginalfft = np.multiply(hpmaskfft,cm)
    print(TheOriginalfft[256,256])
    return TheOriginalfft
    


def power_spectra_window(im,indices,apix):
    import numpy as np
    #indices: startx endx, starty, endy
    #1. FFT
    fft_im = fftim(square(im[indices[0]:indices[1],indices[2]:indices[3]]))
    #2. crop FFT to 3 Angstrom
    #ndx = np.dot( fft_im.shape, apix/3 )
    #c = np.dot( fft_im.shape , 0.5 )
    #s = c - ndx
    #e = c + ndx
    #s = np.ndarray.astype(s,int)
    #e = np.ndarray.astype(e,int)
    #fft_im = fft_im[ (s[0]):(e[0]), (s[1]):(e[1])]
    #3. take absolute value of FFT
    fft_im = (np.abs(fft_im))
    sz = fft_im.shape
    center = np.dot(sz,0.5)
    #6b. high-pass filter the fft
    sigma = 100000
    hpmask = gauss2d(sz[0],sigma,center)
    hpmaskfft = np.multiply(abs(fft_im),hpmask)
    #7. mask out central peak
    #gmask=np.abs(np.add( np.dot( gauss2d(sz[0],10,center), -1 ), 1)) 
    cm=np.add( np.dot( circmask(sz[0],10,center), -1), 1)
    TheOriginalfft = np.multiply(hpmaskfft,cm)
    return TheOriginalfft



def lowpass_filt(m):
    #compute the indices in reciprocal space
    #mm = square(m)
    mm = m
    sz = mm.shape
    testmask = gauss2d(sz[0],5000,[sz[0]/2,sz[1]/2])
    #ndx = np.ndarray.astype(np.dot( mm.shape, angpix/resolution ), int)
    fftm = np.fft.fftshift(np.fft.fft2(mm))
    maskedfft = np.multiply(fftm,testmask)
    
    return np.real(np.fft.ifft2(np.fft.ifftshift(maskedfft)))




def add_scale_bar(im,apix,len):
    #0.75 of x & y
    numpix = len / apix
    sz = im.shape
    xs = (int)(sz[0]*0.75)
    print(xs)
    xe = (int)(sz[0]*0.75+numpix)
    print(xe)
    ys = (int)(sz[1]*0.75)
    print(ys)
    ye = (int)(sz[1]*0.75+100)
    print(ye)
    xrange =  np.arange( xs , xe )
    yrange =  np.arange( ys , ye )
    #print(xrange)
    #np.put( im, np.meshgrid(xrange,yrange), 1)
    im[xrange,ys] = 0
    return im






#this is for reading cryosparc .cs files into dataframe
#requires numba, pyfftw
def read_cs(csfile):
    from pyem import metadata
    from pyem import star
    cs = np.load(csfile)
    df = metadata.parse_cryosparc_2_cs(cs)
    return df

def radial_profile(data):
    center = np.dot(data.shape,0.5)
    #from https://stackoverflow.com/questions/21242011/most-efficient-way-to-calculate-radial-profile
    y,x = np.indices((data.shape)) # first determine radii of all pixels
    r = np.sqrt((x-center[0])**2+(y-center[1])**2)
    ind = np.argsort(r.flat) # get sorted indices
    print(len(ind))
    sr = r.flat[ind] # sorted radii
    sim = data.flat[ind] # image values sorted by radii
    ri = sr.astype(np.int32) # integer part of radii (bin size = 1)
    # determining distance between changes
    deltar = ri[1:] - ri[:-1] # assume all radii represented
    rind = np.where(deltar)[0] # location of changed radius
    nr = rind[1:] - rind[:-1] # number in radius bin
    csim = np.cumsum(sim, dtype=np.float64) # cumulative sum to figure out sums for each radii bin
    tbin = csim[rind[1:]] - csim[rind[:-1]] # sum for image values in radius bins
    radialprofile = tbin/nr # the answer
    return radialprofile
