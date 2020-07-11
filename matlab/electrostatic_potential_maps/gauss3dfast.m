rfunction[outp] = gauss3d(Mi,Mj,Mk,center,sigma)

% tic
% 
% shft = exp(- (2*pi)*i*( Mi.*center(1) + Mj.*center(2) + Mk.*center(3) ));
% ampl = 1/(sqrt((2*pi)^3))*1/(sigma^3);
% 
% toc
% tic
% 
% outp = ifftn( ifftshift( ampl * exp( -pi^2*(Mi.^2+Mj.^2+Mk.^2)./(2*sigma^2)) .*shft  ));
% 
% toc

ampl = 1/(sqrt((2*pi)^3))*1/(sigma^3);
outp = ifftn( ifftshift( ampl * exp( -pi^2*(Mi.^2+Mj.^2+Mk.^2)./(2*sigma^2) - ( (2*pi)*i*( Mi.*center(1) + Mj.*center(2) + Mk.*center(3) ))  )));

end