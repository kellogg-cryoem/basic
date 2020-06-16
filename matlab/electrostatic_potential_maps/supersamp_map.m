function[supermap]=supersamp_map(m,factor)
    fftm = fftshift( fftn ( m ));
    supfft = padarray(fftm,size(m)*(factor/2),0,'both');
    size(fftm)
    size(supfft)
    supermap = real(ifftn(ifftshift(supfft)));
end