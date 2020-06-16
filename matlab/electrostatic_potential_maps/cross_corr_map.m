function[CC] = cross_corr_map(map1,map2)
    %not sure if these equations are correct??
    convfft = fftshift(fftn(map1)) .* fftshift(fftn(map2));
    CC = sum(sum(sum(convfft,1),2),3);
end