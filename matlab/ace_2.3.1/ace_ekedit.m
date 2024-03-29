function   ctfparams = ace(mrc,stig,medium,scope_params); 
%
% DESCRIPTION: 
%
% Returns the parameters of the Contrast Transfer Function. This function 
% is at the heart of Automated CTF Estimation and is called by both 
% leginon_ace_gui and acedemo. 
%
% USAGE  : 
%        ctfparams = ace(filename,outfile,display,stig,medium,dforig,tempdir);
%
%        filename : Name of the input mrc file.
%        outfile  : The name of output file. 
%        display  : Boolean input to control the graphical display.
%        stig     : Boolean switch to turn on astigmatism estiamtion.
%        medium   :  'carbon' or 'ice'.
%        dforig   : The defocus set by the microscope. It is used only 
%                   when the edge detection fails. It defaults to -2um. 
%        tempdir  : Optional string argument to specify directory for
%                   temporary files. 
%              
%        ctfparams: Vector of ctf parameters. First two element give 
%                   the defoci, the third value is Amplitude contrast, 
%                   the forth value is the angle of astigmatism, elements 
%                   5-8 are the noise parameters and 9-12 are the envelope 
%                   parameters.
%       
% 
% See also leginon_ace_gui and acedemo. 
%
% Copyright 2004-2005 Satya P. Mallick. 

if nargin<4
    stig = 1;
    medium = 'carbon';
    dforig = 2e-6;
    tempdir = './';
elseif nargin<5
    medium = 'carbon';
    dforig = 2e-6;
    tempdir = './';
elseif nargin<6
    medium = 'carbon';
    dforig = 2e-6;
    tempdir = './';
elseif nargin<7
    tempdir = './';
end

warning off all
V = scope_params(1);
Cs = scope_params(2);
Ca = scope_params(3);

V = V*1e3;
Cs = Cs*1e-3;
Ca = Ca; %ang per pix scaling
lambda = getlambda(V);

% End Load Microscope Parameters

%START: Initialization

conf = 0;

%End: Initialization
file = mrc;
filesz = size(file);
file = file - mean(file(:));

if(strcmp(medium,'ice'))
    trial = 3;
    confth = 0.1;
    pfice = 0.300;
    powerfactor = pfice;
        
    
elseif(strcmp(medium,'carbon'))
    confth = 0.1;
    pfcarbon = 0.9;
    powerfactor = pfcarbon;
   
    
else
    fprintf('Medium should be ice or carbon\n');
    return;
end
  
  %START: Get 1D profile 
  
  [val, pg, imfftabs, rat, ang] = getprofile(file,stig,medium,tempdir);

  %END: Get 1D profile 
  
  [imheight, imwidth] = size(imfftabs); 
  freqfactor = (imwidth*Ca)/1e10; 
  
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   pg1 = smooth(pg,round(5*imwidth/512));
   dpg1 = diff(pg1); 
   ddpg1 = diff(dpg1); 
 
  
  %START:  Determining the lower cuttoff 
  
  load(strcat(tempdir,'k1')); 
  load(strcat(tempdir,'k2')); 
  x = round(max(sqrt([1/k1 1/k2])));
  
  
  lambda = getlambda(V); 
  
  xtemp = sqrt((-abs(dforig)*lambda+sqrt(dforig^2*lambda^2+2*Cs*lambda^3))...
      /(Cs*lambda^3)); 
  xtemp = round(xtemp*freqfactor); 
  
  if(abs(xtemp-x)>imwidth/20)
    x = xtemp-round(imwidth/50); 
  end 
  x1 = x;
  
   if(~isreal(x) | x>0.35*imwidth)
   
    fprintf(outid,'%s %s\n', filename,...
	'Could not process micrograph:Low SNR'); 
    fprintf('%s %s\n', filename,...
	'Could not process micrograph:Low SNR '); 
    return;
   end  
  
  
  
  %END:  Determining the lower cuttoff 
   
  %START: Removing lower frequencies. 
  
  pgram = pg1(x:end); 
  pgram1 = pg1(x1:end); 
  
  %END: Removing lower frequencies. 
  
  %START: Calculation of the noise function
    
  s = [x:length(pg1)]';
  pixeloff = 1; 
  s = s - pixeloff; 
   
    
  snoise = [x1:length(pg1)]';
  snoise = snoise - pixeloff; 
    
    
  A1 = [ones(length(snoise),1) sqrt(snoise) snoise snoise.^2];
  
  b1 = log(pgram1); 
  %minimize objective function for fit to ctf
   options1 = optimset('TolFun',1e-4,'MaxFunEvals',...
      30,'MaxIter',1000,'Display','off'); 
  [optparams1,fun1,eflag1,out1] = fmincon(@objnoise,[0 0 0 0]', A1,b1,...
      [],[],[],[1e6 0 0 0],[],options1,A1,b1);
  
  if(eflag1<0) 
    fprintf(outid,'%s %s\n', filename,...
	'Could not process micrograph,:Unreliable noise fit'); 
    fprintf('%s %s\n', filename,...
	'Could not process micrograph:Unreliable noise fit'); 
    fclose(outid); 
    return
  end 
  A1 = [ones(length(s),1) sqrt(s) s s.^2];
  b_calc1 = 2*A1*optparams1;
    
  noiseparams = 2*[optparams1(1) optparams1(2)*sqrt(freqfactor) ...
	optparams1(3)*(freqfactor) optparams1(4)*(freqfactor).^2]; 
    
  %END: Calculation of the noise function
    
  %START: Calculation of the ctf + env  function 
  
  ctfenv = pgram.^2 - exp(b_calc1); 
  
  %END: Calculation of the ctf + env  function
  
  %START: Calculation of higher cutoff frequency
  
  ctfenvcum = cumsum(ctfenv(x1-x+1:end)); 
  upcut = find(ctfenvcum> powerfactor*ctfenvcum(end)); 
  
  
  if(strcmp(medium,'ice')) 
    upcut=60;  
  end 
  
  %END: Calculation of higher cutoff frequency
    
  ctfenv = ctfenv(1:upcut); 
  s = s(1: upcut);
  pg = pg(x:end); 
  pgram = pg1(x:x+upcut-1); 
  A1 = [ones(length(s),1) sqrt(s) s s.^2];
  b_calc1 = 2*A1*optparams1;
  
  %START: Calculation of the envelope  function 
    
  A2 = [ones(upcut,1) sqrt(s(1:upcut)) s(1:upcut) s(1:upcut).^2 ];
  b2 = log(ctfenv(1:upcut) -min(ctfenv(1:upcut))+1); 
  options2 = optimset('TolFun',1e-4,'MaxFunEvals',30,...
      'MaxIter',1000,'Display','off'); 
  [optparams2,fun2,eflag2,out2] = fmincon(@objenv,[max(b2)   0  0 0]',...
      -A2, -b2,[],[],[],[1e6 0 0 0],[],options2,A2,b2);
   
  envparams = [optparams2(1) optparams2(2)*sqrt(freqfactor) ...
	optparams2(3)*(freqfactor) optparams2(4)*(freqfactor).^2]; 
    
  if(eflag2<0) 
    fprintf(outid,'%s %s\n', filename,...
	'Could not process micrograph:Unreliable envelope fit'); 
    fprintf('%s %s\n', filename,...
	'Could not process micrograph :Unreliable envelope fit'); 
    fclose(outid);
    return
  end 
  A2 = [ones(length(s),1) sqrt(s) s s.^2 ];
  b_calc2 = A2*optparams2; 
  env = exp(b_calc2);     
  % END: Calculation of the envelope  function 
    
  % START: Calculation of the ctf function 
  ctf = ctfenv./env; 
  
  %Dont user fir filters. Not every user will have the signal processing
  %toolbox. 
  %h = fir1(round(8*imwidth/512),[0.0001 0.1]);
  %ctffilt = conv(ctf,h);
  %ctffilt = ctffilt(1+(length(h)-1)/2:end-(length(h)-1)/2);
 
  ctffilt = smooth(ctf,round(5*imwidth/512));
  
   % END: Calculation of the ctf function     
 
  % START:Calculating the location of zeros of ctf 
  ind= []; 
  ctfsign = sign(diff(ctffilt)); 
  for i=1:length(ctfsign)-1;
    if(ctfsign(i)==-1 & ctfsign(i+1)== 1)
      ind = [ind' i+1]';
      end 
    end 
    if(length(ind)==0) 
      ind = 1; 
    end 
    % End:Calculating the location zeros of ctf 
      
      
      % Start: Robust initial estimate of defocus
      
      soffset = x-1; 
      snew = (s)/freqfactor; 
      sind = ind + soffset - pixeloff; 
      sind = sind/freqfactor; 
      
      flag = 1; 
      n=1 ; 
      nzero = []; 
      
      i = [1:length(ind)]'; 
      imat = repmat(i,1,length(ind)); 
      sindmat = repmat(sind,1,length(ind)); 
      nzero = []; 
      nzero  = (imat -Cs*lambda^3*sindmat.^4/2)./(lambda*sindmat.^2);
     
      err = []; 
      for i=1:length(ind)
	for j=1:length(ind) 
	  gmma = getgamma(snew,nzero(i,j),Cs,lambda); 
	  ctf_temp = getctf(gmma,0.0); 
	  err(i,j) =  norm(ctf-ctf_temp,2); 
	end 
      end
      
      if(sum(isnan(err(:)))==0 &&  sum(isinf(err(:)))==0 && length(err(:))>0)
	[score_final maxscoreind] = min(err(:)); 
	if(strcmp(medium,'carbon'))
	  zinit = nzero(maxscoreind); % initial estimate of defocus 
	else 
	  zinit = nzero(1); 
	end 
	Ainit = 0.0; % Setting the initial Amplitude contrast to 0. 
	gmmainit = getgamma(snew,zinit,Cs,lambda); 
	ctfinit = getctf(gmmainit,Ainit); 
	% End : Robust initial estimate of defocus 
	
	% Start: Refined estimate of defocus
    ctffiltnew =ctffilt;
	ctffiltnew(ind)=0; 
	
	defocusoptions = optimset('TolX',1e-10,'MaxIter',30,'MaxFunEvals',1000,'Display','off');
   
    if strcmp(medium,'carbon')

        [param,fun3,eflag3,out3] = fmincon(@defocusobj,[zinit Ainit],[],[],...
            [],[],[0.0 0.0],[10e-6 0.2],[],...
            defocusoptions,ctffiltnew,snew,V,Cs);
    else
      
      [param,fun3,eflag3,out3] = fmincon(@defocusobj,[zinit Ainit],[],[],...
	  [],[],[zinit- 0.1*1e-6 0.0],[zinit+0.1*1e-6 0.2],[],...
	  defocusoptions,ctffiltnew,snew,V,Cs); 
    end 
	

	if(eflag3<0) 
	  
	  fprintf(outid ,'%s %s\n', filename,...
	      'Could not process micrograph:Unreliable defocus fit'); 
      
	  fprintf('%s %s\n', filename,...
	      'Could not process micrograph:Unreliable defocus fit'); 
	  fclose(outid); 
	  return
  
    end 
	
      	zfinal = param(1); 
	Afinal = param(2); 
	gmmafinal = getgamma(snew,zfinal,Cs,lambda); 
	ctffinal = getctf(gmmafinal,Afinal); 
	% End: Refined estimate of defocus
      
	% End: Estimation of B-factor 
       
	% Start: Calculating the confidence of estimation 
	 meanctffinal = mean(ctffinal); 
	 meanctf_filt = mean(ctffilt); 
	       
	 conf = mean((ctffinal - meanctffinal).*(ctffilt -meanctf_filt)); 
	 conf = conf/(std(ctffinal)*std(ctffilt)); 
     
    end 
     
    % End: Calculating the confidence of estimation 


  slarge =sqrt((-zfinal*lambda+sqrt(zfinal^2*lambda^2+2*Cs*lambda^3))...
      /(Cs*lambda^3)); 
  ssmall = slarge*rat; 
  zsmall = (1-Cs*lambda^3*ssmall^4/2)/(lambda*ssmall^2); 
  defoci = [zfinal zsmall];   Ampconst = Afinal; 
  lowercutoff = (x-pixeloff)/freqfactor;
  uppercutoff = (x+upcut-pixeloff)/freqfactor;

   ssnr = [round(imwidth/20):round(imwidth/4)]'; 
   ssnr = ssnr*1e10/(imwidth*Ca); 
   Asnr = [ones(length(ssnr),1) sqrt(ssnr) ssnr ssnr.^2];
   noise_snr = exp(Asnr*noiseparams'); 
   env_snr   = exp(Asnr*envparams'); 
   snr= sum((env_snr./noise_snr));
  
  ctfparams = [defoci(1), defoci(2) zinit Ampconst -ang*pi/180, noiseparams/2,...
   envparams/2 lowercutoff uppercutoff snr conf];

end