addpath('~/src/Matlab_EM/basic/')
addpath('~/src/Matlab_EM/EMIODist/')
addpath('~/src/Matlab_EM/MRCIO/')

pdb = pdbread('test/TYR-12-B.pdb');
mrc = ReadMRC('test/TYR-12-B.mrc');
%crop before scaling

factor = 4;
supermrc = supersamp_map_2(mrc,factor);
size(mrc)
size(supermrc)

orig_apix = 0.656;
apix = orig_apix/factor;

[mrc,pdb]=crop_map_2(supermrc,pdb,apix);
%make intensity distributions roughly match
mrc(find(mrc < 0))=0;
mrc = norm_mat(mrc);

etbl = read_etable('etable_def.m');

test_std = [0.5:0.5:4];

sumdiffs = zeros(1,8);

for(ii = 1:length(test_std))
    subplot(2,8,ii);
    etbl('C') = [0.5, test_std(ii)];
    etbl('N') = [0.5, test_std(ii)];
    etbl('O') = [0.5, test_std(ii)];
    
    simmap = norm_mat(sim_map_inplace(pdb,mrc,apix,etbl));
    %not needed for now
    %solventmask = sigworth_solvation_mask(pdb,mrc,apix,etbl);
    %simmap = (norm_mat(simmap)+norm_mat(solv)*0.1)*0.5;
    %simmap = simmap + 0.25*solventmask;
    showImage(sum(simmap,3));
    pause(1)
    diff = (sqrt(diff_map(mrc,simmap).^2));
    subplot(2,8,ii+8);
    showImage(sum(diff,3));
    pause(1)
    sumdiffs(ii) = sum(sum(sum(diff)));
    WriteMRC(simmap,apix,strcat('testsimmap-',num2str(test_std(ii)),'.mrc'));
    WriteMRC(mrc,apix,'expmap.mrc');
    WriteMRC(diff,apix,strcat('diffmap-',num2str(test_std(ii)),'.mrc'));
end
