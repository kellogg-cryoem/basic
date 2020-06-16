addpath('~/src/Matlab_EM/MRCIO');
addpath('~/src/Matlab_EM/EMIODist');
addpath('~/src/Matlab_EM/basic/');
%tyr = ReadMRC('tyr29A.mrc');
%tt = tyr(135:150,70:100,190:220);

tyr = ReadMRC('tyr-cropped2.mrc');
contour3((squeeze(sum(tyr,3))))
%need to supersample matrix to see contours more smoothly
tyr_rot = interp3(imrotate3(tyr,100,[1,0,0]));
for(i = 22:46)
subplot(5,5,j)
j=j+1;
contour(squeeze(tyr_rot(i,:,:)),25)
end

simmap = ReadMRC('molmap_1.5ang.mrc');
%use tricubic interpolation for best accuracy
simmap_rot = interp3(imrotate3(simmap,100,[1,0,0]));
%read in atoms and show as points on plot

%next need to write code that orients sidechain with respect to reference?

%just crop out sidechain only with chimera
%align sidechain (with map) to reference sidechain
%resample onto original box (or better yet just save transform and original
%map)
%then crop boxes to only non-zero values

