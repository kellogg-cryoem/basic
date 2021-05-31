dirname = '/ctf/04mar04series/';
list = dir(strcat(dirname,'matfiles/*.mat')); 
for i=1:length(list) 
  load(strcat(dirname,'matfiles/',list(i).name));
  if(ctfparams~=-1)
    if(abs(abs(dforig)-mean(ctfparams(1:2)))<1e-6)
      s = [25:100]'; 
      s = s*1e10/(512*scopeparams(3)); 
      A = [ones(length(s),1) sqrt(s) s s.^2];
      noise = exp(2*A*ctfparams(5:8)'); 
      env   = exp(2*A*ctfparams(9:end)'); 
      snr= sum((env./noise));
      dbctf(list(i).name(1:end-4),dforig,ctfparams,strcat(dirname,'opimages/',list(i).name(1:end-4),'1.jpg'),strcat(dirname,'opimages/',list(i).name(1:end-4),'2.jpg'),snr);
      fprintf(strcat(list(i).name(1:end-4),' written\n'));
    end 
  end   
end
