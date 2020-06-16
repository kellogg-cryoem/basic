function[com] = com_map(mrc)
    ind = find(mrc);
    [i,j,k] = ind2sub(size(mrc),ind);
    com = [mean(i),mean(j),mean(k)];
end