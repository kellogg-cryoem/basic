%change path to point to where your downloaded code lives
addpath('~/src/basic/matlab/basic/')
addpath('~/src/basic/matlab/EMIODist/')
addpath('~/src/basic/matlab/MRCIO/')

%change to match your micrograph of interest and angpix scaling
mrc = ReadMRC('test.mrc');
angpix = 1.33;
%adjust these to crop to a region with lots of filaments
crop = [1 2000];
figure()
showImage(mrc(crop,crop));

[f,r,i] = oneDpowerSpectrum(mrc,angpix);
figure()
%the index where you see a peak (that is not part of the power spectra) is
%your layer-line location
%use the formula: 
% layer_line_spacing = (boxsize * angpix)/ (i - (boxsize/2+1))
%where i is the index where the layer-line peak is observed
plot(r,i);


