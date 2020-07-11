function[outmap] = sim_map(atomic_model,ii,jj,kk,apix,etable)
    outmap = zeros(size(ii));

    natoms = length(atomic_model.Model.Atom);

    tic

     for(ndx = 1:natoms )
         atomName = atomic_model.Model.Atom(ndx).AtomName;
         atomName = atomName(1);
        
         if( isKey(etable,atomName) )
             val = etable(atomName);
             ii_amp = val(1);
             ii_sigma = val(2);
             %[ii_amp, ii_sigma] = lookup_factors(atomName,etable);
             %pdb to image coordinates
             %pdb / apix 
             ii_center = [atomic_model.Model.Atom(ndx).X,...
                          atomic_model.Model.Atom(ndx).Y,...
                          atomic_model.Model.Atom(ndx).Z];
             %in pixels
             %recenter scaled map.. is there a better way to do this?
             outmap = outmap + ...
                        ii_amp * gauss3dfast(ii,jj,kk,(ii_center'),1/ii_sigma);
         end
     end
%%    
    toc

end