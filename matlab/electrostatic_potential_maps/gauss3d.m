function[outp] = gauss3d(inp,center,sigma)
    exponent = sum((inp-center).^2)./(2*sigma^2);
    ampl = 1/ (sigma * sqrt(2*pi));
    outp = ampl*exp(-exponent);
end