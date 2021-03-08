#!/usr/bin/env python
# coding: utf-8

# This piece of code calculates a CTF for a user-defined particle image, and is oriented for use in the Kellogg lab FSC calculator script. This code refers to the code used in the MRC lab CTF calculator script developed by Takanori Nakane https://3dem.github.io/relion/ctf.html
import math
import numpy as np

def conv_electron_wv(U): #Use this function to convert wavelength from V to Angstrom
    h = float(6.626*math.pow(10,-34)) # Planck constant [Js]
    e = float(1.602*math.pow(10,-19)) # electron charge [C]
    c = float(3.000*math.pow(10,8)) # speed of light [m/s]
    m0 = float(9.109*math.pow(10,-31)) # electron rest mass [kg]

    return h / math.sqrt(2 * m0 * e * U) / math.sqrt(1 + e * U / (2 * m0 * c * c)) * 1E10

def get_ctf(apix,boxsize,defocusU,defocusV,astigangle,cs,energy,ampcont,phase_shift):
    apix = apix
    boxsize = boxsize
    defocus = (defocusU+defocusV)/2
    astigmatism = (defocusU - defocusV)/2
    astig_angle = (astigangle/180.0)*math.pi
    cs = cs #spherical abberation
    ht =  energy * 1e3 #energy
    ac = ampcont #amplitude contrast
    phase_shift = (phase_shift/180.0)*math.pi
    
    power_spectrum = False
    img_data = np.zeros([boxsize,boxsize])
    
    wavelength = conv_electron_wv(ht)
    wavelength3 = wavelength*wavelength*wavelength
    pc = math.sqrt(1-(ac*ac))
    K1 = (math.pi/2)*cs*wavelength3
    K2 = math.pi * wavelength
    
    idx = 0
    
    for i in range(0,boxsize):
        sy = (i-(boxsize/2.0))/boxsize/apix
        sy2 = sy*sy
        for j in range(0,boxsize):
            sx = (j-(boxsize/2))/boxsize/apix
            angle = math.atan2(sy,sx)-astig_angle
            local_defocus = defocus +astigmatism*math.cos(2*angle)
            
            s2 = sy2 + (sx**2)
            gamma = K1*s2*s2-K2*s2*local_defocus - phase_shift
            ctf = -pc * math.sin(gamma) + ac * math.cos(gamma)
            
            #color = 0
            #if (power_spectrum == True):
           #     ctf = ctf**2
            #    color = ctf*256
           # else:
           #     color = (ctf+1)*128
            
           # if(color < 0):
           #     color = 0
           # if(color > 255):
          #      color = 255
                
            
            img_data[i,j] = ctf
    return img_data
