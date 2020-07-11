addpath('~/src/Matlab_EM/basic/')
addpath('~/src/Matlab_EM/EMIODist/')
addpath('~/src/Matlab_EM/MRCIO/')




%EMD10101 0.656
%EMD9914  0.7858
%scheres_apoF 0.4

%emdids = {'EMD10101', 'EMD9914', 'scheres_apoF'};
%vals =  [0.656, 0.7858, 0.4]./4;

%angpix_ii = containers.Map(emdids,vals);
%niter = 20;
angpix_ii = 0.656/4;

etbl = read_etable('etable_def.m');

test_std = [0.05:0.1:1.5];
test_amp = [0.05:0.1:1.5];

maps = {};
pdbs = {};

txt = textread('lst','%s');
apix = 0.656/4;
ii = 1;

for(ii = 1:length(txt))
    itm_ii = txt{ii};
    %delim = strfind(itm_ii,'/');
    %dir = itm_ii(1:delim(end));
    %key = strrep( strrep( itm_ii(1:delim(2)), './' , ''), '/', '');
    key = itm_ii;
    prefix = itm_ii(delim(end)+1:end);
    
    pdb = pdbread(strcat(itm_ii,'.pdb'));
    mrc = ReadMRC(strcat(itm_ii,'.mrc'));
    maps{ii} = mrc;
    pdbs{ii} = pdb;
    %apix(ii) = angpix_ii(key);
    
end

kys = keys(etbl);

for(ii = 0:8)
    %pick a random atom to sample
    %atom_ii = kys{randi(length(kys))};
    atom_ii = kys{mod(ii,4)+1};
    objfunc = zeros(length(test_std),length(test_amp));
    for(jj = 1:length(test_std))
        
        for(kk = 1:length(test_amp))
           
           etbl(atom_ii) = [test_amp(kk), test_std(jj)];
           objfunc_ii_kk = 0;
           
           for(ll = 1:length(pdbs))
               
               simmap = norm_mat(sim_map_inplace(pdbs{ll},maps{ll},apix(ll),etbl));
               diff = (sqrt(diff_map(maps{ll},simmap).^2));
               objfunc_ii_kk = objfunc_ii_kk+sum(sum(sum(diff)));
%               WriteMRC(simmap,apix,strcat(txt{ll},'-',num2str(test_std(ii)),'simmap.mrc'));
%               WriteMRC(diff,apix,strcat(txt{ll},'-',num2str(test_std(ii)),'diffmap.mrc'));

           end
           
           objfunc(jj,kk) = objfunc_ii_kk;
           
        end
    end
    ii
    save(strcat('output/objfunc-',atom_ii,'-',num2str(ii),'.mat'),'objfunc');
    save(strcat('output/etbl-',atom_ii,'-',num2str(ii),'.mat'),'etbl');
    [rowj,colk] = find(objfunc == min(min(objfunc)));
    etbl(atom_ii) = [test_amp(colk) test_std(rowj)];
    ii
    atom_ii
    min(min(objfunc))
    test_amp(colk)
    test_std(rowj)
end

