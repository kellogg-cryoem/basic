addpath('~/src/Matlab_EM/basic/')
addpath('~/src/Matlab_EM/EMIODist/')
addpath('~/src/Matlab_EM/MRCIO/')

pdb = pdbread('test/TYR-12-B.pdb');
mrc = ReadMRC('test/TYR-12-B.mrc');
%crop before scaling

factor = 2;
supermrc = supersamp_map(mrc,factor);
size(mrc)
size(supermrc)
scale = size(supermrc)./size(mrc);
orig_apix = 0.656;
apix = orig_apix/scale(1);

[mrc,pdb]=crop_map(supermrc,pdb);

etbl = read_etable('etable_def.m');

test_amp = [0.5:0.5:4];

sumdiffs = zeros(1,8);

for(ii = 1:length(test_amp))
    subplot(2,8,ii);
    etbl('C') = [test_amp(ii), 1];
    etbl('N') = [test_amp(ii), 1];
    etbl('O') = [test_amp(ii), 1];
    
    simmap = sim_map_inplace(pdb,mrc,apix,etbl);
    simmap = align_maps(simmap,mrc);
    showImage(sum(simmap,3));
    pause(1)
    diff = (sqrt(diff_map(mrc,simmap).^2));
    subplot(2,8,ii+8);
    showImage(sum(diff,3));
    pause(1)
    sumdiffs(ii) = sum(sum(sum(diff)));
    %WriteMRC(simmap,apix,'testsimmap.mrc');
%WriteMRC(mrc,apix,'expmap.mrc');
end
