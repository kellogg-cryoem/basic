function[supermap]=supersamp_map(m,factor)
%input matrix has to have even dimensions and factor has to be even 
    fftm = fftshift( fftn ( m ));
    padding = ( size(m)*factor - size(m) ) ./2;
    supfft = padarray(fftm,padding,0,'both');
    supermap = real(ifftn(ifftshift(supfft)));
end