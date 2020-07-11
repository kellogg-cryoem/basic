pdb = pdbread('chainA.pdb');
mrc = ReadMRC('chainA.mrc');
ticks = 1:length(mrc);
[X,Y,Z] = meshgrid(ticks,ticks,ticks);
len1dvec = size(X,1)*size(X,2)*size(X,3); 
meshvec = [reshape(X,1,len1dvec);...
           reshape(Y,1,len1dvec);...
           reshape(Z,1,len1dvec)];

%tt=gauss3d(meshvec,[100,200,300]',2);
etbl = read_etable('etable_def.m');
apix = 0.656;

realsim_pdb = sim_map_inplace(pdb,mrc,apix,etbl);

%freq domain spacing is off
box = 349;
ij = (fftfreq(box,1./(2*apix) ));

ii=zeros(box,box,box);
jj=zeros(box,box,box);
kk=zeros(box,box,box);

for(ind = 1:box)
    for(j = 1:box)
        for(k = 1:box)
            ii(ind,j,k) = ij(ind);
            jj(ind,j,k) = ij(j);
            kk(ind,j,k) = ij(k);
        end
    end
end

ii = fftshift(ii);
jj = fftshift(jj);
kk = fftshift(kk);

tic

fouriersim_pdb = sim_map_inplace_fast(pdb,ii,jj,kk,apix,etbl);
%Mm = gauss3dfast(ii,jj,kk,[100,200,300],1/4);

toc


subplot(1,2,1)
showImage(squeeze(sum( realsim_pdb, 1)));
subplot(1,2,2)
showImage(squeeze(sum( real(fouriersim_pdb), 1)));