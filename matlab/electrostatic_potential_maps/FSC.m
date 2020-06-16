function[fsc] = FSC(map1,map2)
    fft1=fftshift(fftn(map1));
    fft2=fftshift(fftn(map2));
    %complex conjugate
    ccfft2=conj(fft2);
    %https://en.wikipedia.org/wiki/Fourier_shell_correlation
    
    num = dot(fft1,ccfft2);
    sumsqfft1 = abs(fft1).^2;
    sumsqfft2 = abs(fft2).^2;
    
    
end