function[mc,pc] = crop_map(m,p,apix)

%row -> Y
%col -> X

%pdb coordinates in angstrom. image in voxels
% pdb ./ apix -> image
% image * apix -> pdb

limits = zeros(3,2);
%find non-zero indices of refmap
ind = find(m > 1e-3);
[i,j,k] = ind2sub(size(m), ind);
limits = [min(i), max(i);...
          min(j), max(j);...
          min(k), max(k)];

newdim = max([limits(1,2)-limits(1,1)+1,...
             limits(2,2)-limits(2,1)+1,...
             limits(3,2)-limits(3,1)+1]);
           
if(mod(newdim,2) == 1)
    newdim = newdim + 1;
end

centers = round(mean(limits,2));
lb = centers - newdim./2;
ub = centers + newdim./2;
mc = m(lb(1):ub(1),...
          lb(2):ub(2),...
          lb(3):ub(3));

pdbX = ([p.Model.Atom(:).X]) - ((lb(1)-2)*apix); 
pdbY = ([p.Model.Atom(:).Y]) - ((lb(2)-2)*apix);
pdbZ = ([p.Model.Atom(:).Z]) - ((lb(3)-2)*apix);

pc = p;
nat = length( p.Model.Atom(:) );
for(i = 1:nat)
    pc.Model.Atom(i).X = pdbX(i);
    pc.Model.Atom(i).Y = pdbY(i);
    pc.Model.Atom(i).Z = pdbZ(i);
end
end