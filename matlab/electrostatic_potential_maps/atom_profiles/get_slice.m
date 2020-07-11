function[slice] = get_slice(atom1,atom2,view_dir,map,apix,thickness)
   bnd = atom2-atom1;
   ortho_vec = cross(bnd,view_dir);
   indices = ceil( thickness / apix );
   slice = zeros(size(obliqueslice(map,(atom1)./apix,ortho_vec./apix,'OutputSize','full')));
   for(ii = -(indices/2):1:indices/2 )
       slice = slice + obliqueslice(map,(atom1+ii)./apix,ortho_vec./apix,'OutputSize','full');
   end
end