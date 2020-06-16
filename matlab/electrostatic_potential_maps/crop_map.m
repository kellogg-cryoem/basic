function[mc,pc] = crop_map(m,p)

    limits = zeros(3,2);
    %find non-zero indices of refmap
    ind = find(m > 1e-3);
    [i,j,k] = ind2sub(size(m), ind);
    limits = [min(i), max(i);...
              min(j), max(j);...
              min(k), max(k)];

    cropped_sizes = [limits(1,2)-limits(1,1);...
                     limits(2,2)-limits(2,1);...
                     limits(3,2)-limits(3,1)];
    %limits(1,2) -> 1
    %limits(2,2) -> 1
    %limits(3,2) -> 1
    %so subtract limits from pdb coords
                 
    centers = limits(:,1) + (cropped_sizes .*0.5);
    centers(find(mod(centers,2))) = centers(find(mod(centers,2)))+1;

    bsize = max(cropped_sizes);
    if(mod(max(cropped_sizes),2) == 1)
        bsize = max(cropped_sizes) + 1;
    end
    bsize = bsize * 2;
     ll = [centers-bsize/2,...
           centers+bsize/2];
    mc = m(ll(1,1):ll(1,2),ll(2,1):ll(2,2),ll(3,1):ll(3,2));
    nat = length( p.Model.Atom(:) );
    pc = p;
    
    %new center of map needs to be center of pdb
    pdbcoords = [squeeze([p.Model.Atom(:).X]); ...
                 squeeze([p.Model.Atom(:).Y]); ...
                 squeeze([p.Model.Atom(:).Z])];
    pdb_com = mean(pdbcoords');
    ind = find(mc);
    [i,j,k] = ind2sub(size(mc), ind);
    mc_com = [mean(i), mean(j), mean(k)];
    
    for(i = 1:nat)
        pc.Model.Atom(i).X = p.Model.Atom(i).X - (pdb_com(1)-mc_com(1));
        pc.Model.Atom(i).Y = p.Model.Atom(i).Y - (pdb_com(2)-mc_com(2));
        pc.Model.Atom(i).Z = p.Model.Atom(i).Z - (pdb_com(3)-mc_com(3));
    end
end