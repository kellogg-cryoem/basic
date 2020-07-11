function[outp] = gauss3d(inp,center,sigma)
    exponent = sum((inp-center).^2)./(2*sigma^2); %this is the slowest step
    ampl = 1/ (sigma^3 * (2*pi)^1.5);
    outp = ampl*exp(-exponent); 
end