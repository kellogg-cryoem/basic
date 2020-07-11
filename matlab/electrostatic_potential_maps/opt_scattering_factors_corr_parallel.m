addpath('./Matlab_EM/basic/')
addpath('./Matlab_EM/EMIODist/')
addpath('./Matlab_EM/MRCIO/')
addpath('./src/')



%EMD10101 0.656
%EMD9914  0.7858
%scheres_apoF 0.4

emdids = {'EMD10101', 'EMD9914', 'scheres_apoF'};
vals =  [0.656, 0.7858, 0.4];

angpix_ii = containers.Map(emdids,vals);

etbl = read_etable('etable_def.m');

test_std = [0.05:0.5:1.5];
test_amp = [0.05:0.5:1.5];

maps = {};
pdbs = {};

txt = textread('lst','%s');
apix = zeros(1,length(txt));
ii = 1;

for(ii = 1:length(txt))
    itm_ii = txt{ii};
    delim = strfind(itm_ii,'/');
    dir = itm_ii(1:delim);
    key = strrep( strrep( itm_ii(1:delim-1), './' , ''), '/', '')
    prefix = itm_ii(delim+1:end);
    
    pdb = pdbread(strcat(itm_ii,'.pdb'));
    mrc = ReadMRC(strcat(itm_ii,'.mrc'));
    maps{ii} = mrc;
    pdbs{ii} = pdb;
    apix(ii) = angpix_ii(key)
    
end

kys = keys(etbl);

%create parfor loop, reqs:
%   - n is integer, increasing
%   - no nesting

%enumerate parameters,
numkeys = length(kys);
pars = meshgrid(test_std,test_amp,1:numkeys);

iter = sub2ind(size(pars),length(test_std),length(test_amp),numkeys);
iter
objfuncs = {};
for( i = 1:numkeys )
   objfuncs{i} = zeros(length(test_std),length(test_amp)); 
end

keyset = {'C','N','O','H'};

M = 3; %num workers
gpu_ids = [0 1 2];

for( ii = 1:iter)
    ii
    [ jj,kk,ky_i] = ind2sub(size(pars),ii);
    %spmd
    %   gpuDevice( 1+ mod( labindex - 1, gpuDeviceCount ) )
    %end
    
    %pick a random atom to sample
    vals = {[2 3],[2 3],[2 3],[0.25 2]};  
    
    vals{ky_i} = [test_amp(kk), test_std(jj)];
    %etbl = containers.Map(keyset,vals);
    %   why is struct a problem with parfor?
    %    etbl(atom_ii) = [test_amp(kk), test_std(jj)];
    objfunc_ii_kk = 0;
           
    for(ll = 1:length(pdbs))
       simmap = norm_mat(real(sim_map_inplace_fastgpu(pdbs{ll},maps{ll},apix(ll),etbl)));
       tt = (corrcoef(maps{ll},simmap));
       diff = tt(1,2);
       objfunc_ii_kk = objfunc_ii_kk+diff;
               WriteMRC(simmap,apix,strcat(txt{ll},'-',num2str(test_std(jj)),'-',num2str(test_amp(kk)),'simmap.mrc'));

    end
    oo = objfuncs{ky_i};       
    oo(jj,kk) = objfunc_ii_kk;
    objfuncs{ky_i} = oo;

    save(strcat('output/objfunc-',ky_i,'-',num2str(ii),'.mat'),'objfunc_ii_kk');
    save(strcat('output/etbl-',ky_i,'-',num2str(ii),'.mat'),'etbl');
    [rowj,colk] = find(objfunc_ii_kk == max(max(objfunc_ii_kk)));
    %etbl(ky_i) = [test_amp(colk) test_std(rowj)];
end

