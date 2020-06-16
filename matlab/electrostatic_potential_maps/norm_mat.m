function[nmat] = norm_mat(mm)
    min_mm = min(min(min(mm)));
    max_mm = max(max(max(mm)));
    nmat = ( mm - min_mm )./ (max_mm - min_mm);
end