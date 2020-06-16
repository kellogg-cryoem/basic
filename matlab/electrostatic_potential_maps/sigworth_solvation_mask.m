%take code from sigworth's SolventAndProteinDensity in order to compute
%solvation 

function[s] = sigworth_solvation_mask(atomic_model,inmap,apix,etable)
    %loop over atoms, compute distance to nearest atom.
    natoms = length(atomic_model.Model.Atom);
    coords = [squeeze([atomic_model.Model.Atom(:).X]);...
              squeeze([atomic_model.Model.Atom(:).Y]);...
              squeeze([atomic_model.Model.Atom(:).Z])];
    
    mapdim = size(inmap);
    [i,j,k] = meshgrid(1:mapdim(1),1:mapdim(2),1:mapdim(3));
    i = reshape(i,1,size(i,1)*size(i,2)*size(i,3));
    j = reshape(j,1,size(j,1)*size(j,2)*size(j,3));
    k = reshape(k,1,size(k,1)*size(k,2)*size(k,3));
    %voxels x apix = pdb-coordinates
    dmat = zeros(mapdim);
    for(ii = 1:natoms)
        if(ii == 1)
            for(nn = 1:length(i))
             dmat(i(nn),j(nn),k(nn)) = sqrt((i(nn).*apix - coords(1,ii)).^2 + ...
                                (j(nn).*apix - coords(2,ii)).^2 + ...
                                (k(nn).*apix - coords(3,ii)).^2 );
            end
        else
            for(nn = 1:length(i))
             dmat(i(nn),j(nn),k(nn)) = min(dmat(i(nn),j(nn),k(nn)),sqrt((i(nn).*apix - coords(1,ii)).^2 + ...
                                                (j(nn).*apix - coords(2,ii)).^2 + ...
                                                (k(nn).*apix - coords(3,ii)).^2 ));
            end
        end
    end
    WaterDensity = 4.9327; %taken from SolventAndProteinDensity.m
    %r1-r3 is polar radius
    %take from etable
    %r1 = etable('C');
    %r1 = r1(2)
    %r2 = etable('C');
    %r2 = r2(2);
    %r3 = etable('C');
    %r3 = r3(2);
    r1 = 2.5; r2 = 2.5; r3 = 2.5;
    s=WaterDensity*((1+erf(dmat-r1))./2+0.2.*exp(-((dmat-r2)./2.5).^2)-0.15.*exp(-((dmat-r3)./1.5).^2));

end