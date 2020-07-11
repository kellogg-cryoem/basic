function[c] = com(pdb,nums)
    c = zeros(1,3);
    for(i = 1:length(nums))
        c = c + get_atom_coords(pdb,nums(i));
    end
    c = c ./ length(nums);
end