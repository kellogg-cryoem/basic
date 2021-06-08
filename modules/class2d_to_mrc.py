#!/usr/bin/python


def read_cs(csfile):
    from pyem import metadata
    from pyem import star
    import numpy as np
    cs = np.load(csfile)
    df = metadata.parse_cryosparc_2_cs(cs)
    return df


def get_micrographnames(df):
    micrographnames =[]
    for i in df.rlnImageName:
        micrographnames.append(i.split('/')[-1])
    return micrographnames

def image_indices(lst,nm):
    #assumes nm is a string
    ndc = []
    for nn, nn_name in enumerate(lst):
        if( nn_name == nm ):
            ndc.append(nn)
    return ndc

def unique(lst):
    return list(set(lst))

def count_particles(df):
    #return a list of micrographs and number of particles per micrograph
    mrcs = get_micrographnames(df)
    uniq_mrcs = unique(mrcs)
    numptcls = []
    for ii, ii_mrc in enumerate(uniq_mrcs):
        numptcls.append( len(image_indices( mrcs, ii_mrc )) )
    if( len(numptcls) != len(uniq_mrcs) ):
        print("dont know how this happened but there's a mismatch between your list lengths!")
    return[ uniq_mrcs, numptcls ]

def get_particles(df,micrographname):
     mrcs = get_micrographnames(df)
     nn = image_indices(mrcs,micrographname)
     return df.loc[nn]



def lowpass_filt_2darray(nparr):
    import basicEM
    import numpy as np
    mm = basicEM.square(nparr)
    sz = mm.shape
    testmask = basicEM.gauss2d(sz[0],1000,[sz[0]/2,sz[1]/2])
    #ndx = np.ndarray.astype(np.dot( mm.shape, angpix/resolution ), int)
    fftm = np.fft.fftshift(np.fft.fft2(mm))
    maskedfft = np.multiply(fftm,testmask)
    
    return np.real(np.fft.ifft2(np.fft.ifftshift(maskedfft)))


def copy_paste_offset(src,dest,offset):
    import numpy as np
    import cv2 as cv2
    from matplotlib import pyplot as plt
    #take 75% of area
    
    sshape = src.shape
    c = (np.ceil(np.dot(sshape,0.5)))
    halfc = (np.ceil(np.dot(c,0.5)))
    s = 1.2
    src_crop = src[int(c[0]-(s*halfc[0])):int(c[0]+(s*halfc[0])),int(c[1]-(s*halfc[1])):int(c[1]+(s*halfc[1]))]
    tt = lowpass_filt_2darray(src_crop)

    #medianBlur doesn't seem to work
    #plt.imshow(tt)
    img = cv2.medianBlur(np.dot(np.array(tt),255).astype('uint8'),5)

    ndx = np.argwhere(img < 160)
    #print(ndx)
    #test = np.zeros((100,100))
    for i in range(len(ndx)):
        #print(ndx[i,])
        nn = ndx[i,]
        dest[int(offset[0]+nn[0]),int(offset[1]+nn[1])] = src_crop[nn[0],nn[1]]
    return dest

import scipy

def map2d_classes_to_orig(classavg_imgs,all_coords_x,all_coords_y,mrc00,szx,szy,output_name):
    import matplotlib.pyplot as plt
    import numpy as np
    #fig = plt.figure()
    #ax1 = fig.add_subplot(121)
    #ax1.imshow(basicEM.lowpass_filt(m), interpolation='bilinear', cmap='gray')
    mins = np.min(((int)(szx),(int)(szy)))
    s = (mins,mins)
    background = np.zeros(s)
    print(background.shape)
    import skimage as skimage
    import basicEM

    imgdata = classavg_imgs.data

    for ii,ptcl_ii in enumerate(all_coords_x):
        x = np.ceil(all_coords_x[ii])
        y = np.ceil(all_coords_y[ii])
        ptclii = mrc00.iloc[ii]
        #extents: left, right, bottom, top
        xformimg = scipy.ndimage.interpolation.rotate(imgdata[ ptclii.rlnClassNumber-1, ], ptclii.rlnAnglePsi*-1 )
        pshape = xformimg.shape
        xformimg = skimage.transform.resize(xformimg,np.dot(pshape,4))
        pshape = xformimg.shape
        #print(pshape)
        #scipy.ndimage.shift(xformimg,[0,0],background)
        #print(x)
        #print(x+pshape[0])
        #print(y)
        #print(y+pshape[1])
        if( x+pshape[0] < (int)(szx) and y+pshape[1] < int(szx) ):
            background = copy_paste_offset(xformimg,background,[int(y),int(x)])
            #take the middle 50%
            #c = (np.ceil(np.dot(pshape,0.5)))
            #halfc = (np.ceil(np.dot(c,0.5)))
            #s = 1.2
            #xformimg = xformimg[int(c[0]-(s*halfc[0])):int(c[0]+(s*halfc[0])),int(c[1]-(s*halfc[1])):int(c[1]+(s*halfc[1]))]
            #pshape = xformimg.shape
            #background[int(x):int(x+pshape[0]),int(y):int(y+pshape[1])] = xformimg
            #background[int(x):int(x+240),int(y):int(y+240)] = xformimg[25:85,25:85]
        #plt.imshow( scipy.ndimage.interpolation.rotate(imgdata[ ptclii.rlnClassNumber-1, ], ptclii.rlnAnglePsi*-1 ), extent=[x,x+110,y,y+110] )
        #ax1.imshow( scipy.ndimage.interpolation.rotate(imgdata[ ptclii.rlnClassNumber-1, ], ptclii.rlnAnglePsi*-1 ) )
    plt.rcParams['figure.figsize'] = [100, 100]
    plt.imshow(background, interpolation='bilinear', cmap='gray')
    plt.savefig(output_name)
    plt.close()
    
def main():

    import sys
    if( len(sys.argv) == 1 ):
        print("usage: ",sys.argv[0],"<CLASSAVG_MRC> <CRYOSPARC_CS> <PASSTHROUGH_PTCLS> <IMG_SZ_X> <IMG_SZ_Y> <OUTPUT_DIR")
        exit(1)

    sys.path.append('/home/ehk68/src/basic')
    sys.path.append('/home/ehk68/src/pyem-master/')

    import basicEM
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    import mrcfile as mrc
    import numba
    import PIL as PIL

    #do a loop over many images: names = get_micrographnames(classavgs)
    imgs = sys.argv[1]
    cs = sys.argv[2]
    passthrough = sys.argv[3]
    classavg_imgs = mrc.open(imgs)
    img_sz_x = sys.argv[4]
    img_sz_y = sys.argv[5]
    output_dir = sys.argv[6]
    star_df = read_cs(cs)

    #get particle coordinates from passthrough file
    ptcl_coords = np.load(passthrough)
    

    names = get_micrographnames(star_df)
    for (i,val) in enumerate(np.unique(names)):
    	print(val)
    	nn=(image_indices(names,val))
    	all_coords_x = []
    	all_coords_y = []
    	for ii,ith_ndx in enumerate(nn):
       	     all_coords_x.append( (int)(ptcl_coords[ith_ndx][5]*(float)(img_sz_x)) )
             all_coords_y.append( (int)(ptcl_coords[ith_ndx][6]*(float)(img_sz_y)) )
    	mrc00 = get_particles(star_df,val)
    	pngname = val.replace(".mrc",".png") 
    	map2d_classes_to_orig(classavg_imgs,all_coords_x,all_coords_y,mrc00,img_sz_x,img_sz_y,\
                          output_dir+pngname)


if __name__ == "__main__":
   main()
