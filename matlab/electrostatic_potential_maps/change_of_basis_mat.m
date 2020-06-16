ref1pdb = pdbread('ref1.pdb');
ref1pdb = [squeeze([ref1pdb.Model.Atom(:).X]); ...
           squeeze([ref1pdb.Model.Atom(:).Y]); ...
           squeeze([ref1pdb.Model.Atom(:).Z])];

ref2pdb = pdbread('ref2.pdb');
ref2pdb = [squeeze([ref2pdb.Model.Atom(:).X]); ...
           squeeze([ref2pdb.Model.Atom(:).Y]); ...
           squeeze([ref2pdb.Model.Atom(:).Z])];

ref3pdb = pdbread('ref3.pdb');
ref3pdb = [squeeze([ref3pdb.Model.Atom(:).X]); ...
           squeeze([ref3pdb.Model.Atom(:).Y]); ...
           squeeze([ref3pdb.Model.Atom(:).Z])];

ref1mrc = ReadMRC('ref1.mrc');
ind = find(ref1mrc);
[i,j,k] = ind2sub(size(ref1mrc), ind);
ref1mrc = [mean(i), mean(j), mean(k)];

ref2mrc = ReadMRC('ref2.mrc');
ind = find(ref2mrc);
[i,j,k] = ind2sub(size(ref2mrc), ind);
ref2mrc = [mean(i), mean(j), mean(k)];

ref3mrc = ReadMRC('ref3.mrc');
ind = find(ref3mrc);
[i,j,k] = ind2sub(size(ref3mrc), ind);
ref3mrc = [mean(i), mean(j), mean(k)];

basispdb = [ref1pdb, ref2pdb, ref3pdb];
basismrc = [ref1mrc', ref2mrc', ref3mrc'];

% need to change basis of mrc to pdb
%in order to superimpose mrctest to mrcref
[Rmrc2pdb,tmrc2pdb]=rigid_transform_3D(basismrc, basispdb);

[Rpdb2mrc,tpdb2mrc]=rigid_transform_3D(basispdb, basismrc);