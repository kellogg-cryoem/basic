function[outmap] = sim_map(atomic_model,inmap,apix,etable)
    outmap = zeros(size(inmap));
    %make array of gaussian centers, amplitudes, and sigma
    %then as you iterate through each voxel
    %sum output values of gaussians
    if( size(outmap,1) ~= size(outmap,2) || ...
        size(outmap,2) ~= size(outmap,3) || ...
        size(outmap,1) ~= size(outmap,3) )
        'map has to be cube'
        return
    end
    ticks = 1:length(outmap);
    [X,Y,Z] = meshgrid(ticks,ticks,ticks);
    len1dvec = size(X,1)*size(X,2)*size(X,3); 
    meshvec = [reshape(X,1,len1dvec);...
               reshape(Y,1,len1dvec);...
               reshape(Z,1,len1dvec)];
    sum_of_gauss = [];
    natoms = length(atomic_model.Model.Atom);
    
    output1d = zeros(1,length(meshvec));
    natoms
    
    for(ii = 1:natoms )
        atomName = atomic_model.Model.Atom(ii).AtomName;
        atomName = atomName(1);
       
        if( isKey(etable,atomName) )
            val = etable(atomName);
            ii_amp = val(1);
            ii_sigma = val(2);
            %[ii_amp, ii_sigma] = lookup_factors(atomName,etable);
            %pdb to image coordinates
            %pdb / apix 
            ii_center = [atomic_model.Model.Atom(ii).X./apix,...
                         atomic_model.Model.Atom(ii).Y./apix,...
                         atomic_model.Model.Atom(ii).Z./apix];
            %in pixels
            %recenter scaled map.. is there a better way to do this?
            output1d = output1d + ...
                       ii_amp * gauss3d(meshvec,(ii_center'),ii_sigma);
        end
    end
    
    
    for(i = 1:length(meshvec))
        outmap(squeeze(meshvec(1,i)),...
               squeeze(meshvec(2,i)),...
               squeeze(meshvec(3,i))) = output1d(i);
    end
    
end