function[outmap] = sim_map(atomic_model,map,apix,etable)
    outmap = gpuArray(zeros(size(map)));

    natoms = length(atomic_model.Model.Atom);

    box = length(map);
    ij = fftfreq(box,1);
    
    ii = (zeros(box,box,box));
    jj = (zeros(box,box,box));
    kk = (zeros(box,box,box));
   
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
    
    

    ii = gpuArray(ii);
    jj = gpuArray(jj);
    kk = gpuArray(kk);
    

     for(ndx = 1:natoms )
         atomName = atomic_model.Model.Atom(ndx).AtomName;
         atomName = atomName(1);
        
         if( isKey(etable,atomName) )
             val = etable(atomName);
             ii_amp = val(1);
             ii_sigma = val(2);
             s = 1 / ii_sigma;
             %[ii_amp, ii_sigma] = lookup_factors(atomName,etable);
             %pdb to image coordinates
             %pdb / apix 
             ii_center = [atomic_model.Model.Atom(ndx).X./apix,...
                          atomic_model.Model.Atom(ndx).Y./apix,...
                          atomic_model.Model.Atom(ndx).Z./apix];
             %in pixels
             %recenter scaled map.. is there a better way to do this?
             ampl = 1/(sqrt((2*pi)^3))*(1/(s^3));
             
             outmap = outmap + ...
                        ii_amp * ifftn( ifftshift( ampl * exp( -pi^2*(ii.^2+jj.^2+kk.^2)./(2*s^2) - ( (2*pi)*i*( ii.*ii_center(1) + jj.*ii_center(2) + kk.*ii_center(3) ))  )));
             %gauss3dfast(ii,jj,kk,(ii_center'),1/ii_sigma);
         end
     end
%%    
   
 
    
    outmap = gather(outmap);
    

end