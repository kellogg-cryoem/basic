addpath('~/src/Matlab_EM/basic/')
addpath('~/src/Matlab_EM/EMIODist/')
addpath('~/src/Matlab_EM/MRCIO/')
addpath('./src/')

txt = textread('phe.lst','%s');
ii = 1;

%EMD10101 0.656
%EMD9914  0.7858
%scheres_apoF 0.4

keys = {'EMD10101', 'EMD9914', 'scheres_apoF'};
vals =  [0.656, 0.7858, 0.4];

angpix_ii = containers.Map(keys,vals);

for(ii = 1:length(txt))
    itm_ii = txt{ii};
    delim = strfind(itm_ii,'/');
    dir = itm_ii(1:delim(end));
    key = strrep( strrep( itm_ii(1:delim(2)), './' , ''), '/', '');
    prefix = itm_ii(delim(end)+1:end);
    
    pdb = pdbread(strcat(itm_ii,'.pdb'));
    mrc = ReadMRC(strcat(itm_ii,'.mrc'));
    %crop before scaling

    factor = 4;
    supermrc = supersamp_map_2(mrc,factor);
    size(mrc)
    size(supermrc)

    orig_apix = angpix_ii(key);
    %orig_apix = 0.656;
    apix = orig_apix/factor;

    [mrc,pdb]=crop_map_2(supermrc,pdb,apix);
    %make intensity distributions roughly match
    mrc(find(mrc < 0))=0;
    mrc = norm_mat(mrc);

    outmrc = strcat(dir,'/',prefix,'-super.mrc');
    outmrc
    outpdb = strcat(dir,'/',prefix,'-super.pdb');
    
    WriteMRC(mrc,apix,outmrc);
    pdbwrite(outpdb,pdb);
    
    %WriteMRC(simmap,apix,strcat('testsimmap-',num2str(test_std(ii)),'.mrc'));
    %WriteMRC(mrc,apix,'expmap.mrc');
    %WriteMRC(diff,apix,strcat('diffmap-',num2str(test_std(ii)),'.mrc'));
end
