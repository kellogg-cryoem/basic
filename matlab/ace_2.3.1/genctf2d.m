function onlyctf = genctf2d(imwidth,ctfparams,scopeparams) 
%
% DESCRIPTION: 
%     Generates the 2D CTF. 
%
% USAGE: 
%     ctf = genctf2d(imwidth, ctfparams,scopeparams) 
%
%     imwidth     : Width of the image. 
%     ctfparams   : Parameter vector of the CTF parameters. 
%     scopeparams : Parameter vector of microscope parameters. 
%     ctf         : The contrast transfer function. 
%
% Copyright 2004-2005 Satya P. Mallick 

defoci = ctfparams(1:2); 
A = ctfparams(3); %what is A? some coefficient?
ast_ang =  pi/2+ctfparams(4); 
%noisep = ctfparams(5:8); %noise power spectrum?
%envp = ctfparams(9:12);  %envelope function coefficients?

V = scopeparams(1)*1e3; 
Cs = scopeparams(2)*1e-3; 
Ca = scopeparams(3)*1e-10; %ang per pix

if(mod(imwidth,2)==0) %if even sized image
 zro = imwidth/2+0.5; 
else 
  zro = ceil(imwidth/2) 
end 

defocus_mean = (defoci(1)+defoci(2))/2; 
defocus_dev = abs(defoci(1)-defoci(2))/2; 
lambda = getlambda(V); %wavelength of electron at given operating voltage

[i,j] = meshgrid(1:imwidth); 

r = ((i-zro).^2 + (j-zro).^2).^(0.5); %radius
s = r./(imwidth*Ca); % fourier frequency
ang = atan((j-zro)./(i-zro)); 
d = defocus_mean + defocus_dev*cos(2*(ang-ast_ang)); 
gmma = squeeze(getgamma(s,d,Cs,lambda)); 
%env = exp(2*(envp(1)+envp(2).*s.^0.5 + envp(3).*s + envp(4).*s.^2)); 
%noise = exp(2*(noisep(1)+ noisep(2).*s^0.5 + noisep(3).*s + noisep(4).*s.^2));
%ctf = env.*(sqrt(1-A^2)*sin(gmma)+A*cos(gmma)).^2+noise;
ctf = (sqrt(1-A^2)*sin(gmma)+A*cos(gmma)).^2;

%envelope function * CTF + noise

onlyctf = (sqrt(1-A^2)*sin(gmma)+A*cos(gmma));
%[indx1 indy1] = find(abs(gmma-pi)<0.04); 
%save indx1 indx1
%save indy1 indy1




 

