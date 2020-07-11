function[anames] = get_anames(pdb)
    anames = {pdb.Model.Atom.AtomName};
end