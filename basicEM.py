#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import mrcfile as mrc
m = mrc.open('/local1/workdir1/cryosparc/ehk68/P23/J16/motioncorrected/Cascade_grid3_0051_patch_aligned_doseweighted.mrc')


# In[2]:


def square(m):
    #m is a numpy array
    import numpy as np
    import mrcfile 
    mindim = np.min(m.shape)
    return m[0:mindim,0:mindim]


# In[3]:


def normalize(m):
    import numpy as np
    return np.dot(m,1/(sum(sum(m))))


# In[4]:


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


# In[5]:


square_m = square(m.data)
small_m = fft_crop(square_m,0.25)
small_m.shape


# In[6]:


def fftim(m):
    import numpy as np
    return np.fft.fftshift(np.fft.fft2(m))
    


# In[7]:


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

testgauss = np.add( np.dot( gauss2d(512,10,[256,256]), -1 ), 1)
testgauss.shape
plt.imshow(testgauss,vmin=0,vmax=1, cmap='gray' )


# In[8]:


def circmask(sz,radius,center):
    import numpy as np
    ndx = np.linspace(0,sz,sz)
    [X,Y] = np.meshgrid(ndx,ndx)
    mask = (np.square(X-center[0]) + np.square(Y-center[1])) <= radius**2
    mask = 1*mask.astype(float)
    return mask


# In[9]:


cm = circmask(512,10,[256,256])
cm=np.add( np.dot( cm, -1), 1)
plt.imshow(cm)


# In[10]:


cm[0,0]


# In[11]:


testgauss[0,0]


# In[12]:


testmask = gauss2d(512,100000,[256,256])
plt.imshow(testmask)


# In[13]:


def show_power_spectra(im):
    get_ipython().run_line_magic('matplotlib', 'inline')
    from matplotlib import pyplot as plt
    import numpy as np
    plt.rcParams['figure.figsize'] = [25, 25]

    #1. FFT
    fft_im = ((fftim(square(im))))
    #2. crop FFT to 3 Angstrom
    ndx = np.dot( fft_im.shape, 1.28/3 )
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
    


# In[14]:




#minval = np.min(np.min(fft_im))
#maxval = np.max(np.max(fft_im))

#plt.imshow(np.abs(fft_im),vmin=0,vmax=10000000, cmap='gray' )
#plt.show()
#print(minval)
#print(maxval)
#plt.imshow(display_fft(square_m,1.28), interpolation='bilinear', cmap='gray', vmin=minval, vmax=maxval)
#plt.show()
tt = show_power_spectra(square_m)
plt.imshow(abs(tt), interpolation='bilinear',cmap='gray')


# In[15]:


plt.rcParams['figure.figsize'] = [25, 25]
plt.imshow(square(m.data), interpolation='bicubic', cmap='gray', vmin=0, vmax=1)
plt.show()


# In[16]:


tt = np.real(np.fft.ifft2(np.fft.ifftshift(np.fft.fftshift(np.fft.fft2(square(m.data))))))
plt.imshow(tt, interpolation='bicubic', cmap='gray')


# In[39]:


def lowpass_filt(m):
    #compute the indices in reciprocal space
    mm = square(m.data)
    sz = mm.shape
    testmask = gauss2d(sz[0],10000,[sz[0]/2,sz[1]/2])
    #ndx = np.ndarray.astype(np.dot( mm.shape, angpix/resolution ), int)
    fftm = np.fft.fftshift(np.fft.fft2(mm))
    maskedfft = np.multiply(fftm,testmask)
    
    return np.real(np.fft.ifft2(np.fft.ifftshift(maskedfft)))


# In[58]:


np.arange(0,10)


# In[87]:


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


# In[88]:


get_ipython().run_line_magic('matplotlib', 'inline')
from matplotlib import pyplot as plt
plt.rcParams['figure.figsize'] = [25, 25]
print(m.data.shape)
#lpm = fft_crop( lowpass_filt(m), 0.25) #there is a bug in fft_crop
lpm = lowpass_filt(m)
angpix = 1.28
scale_bar_length = 1000 #in angstrom
lpm = add_scale_bar(lpm,angpix,scale_bar_length)
plt.imshow(lpm, interpolation='bilinear', cmap='gray')

