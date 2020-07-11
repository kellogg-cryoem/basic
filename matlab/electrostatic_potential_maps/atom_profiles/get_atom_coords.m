function[coords] = get_atom_coords(pdb,atomNum)
    coords = [pdb.Model.Atom(atomNum).X pdb.Model.Atom(atomNum).Y pdb.Model.Atom(atomNum).Z];
end