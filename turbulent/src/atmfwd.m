function [ifg_atm]=atmfwd(ifglist,dem,refx,refy,stackmethod,atmparams)

nifgs=length(ifglist);
nepoch = max(max(ifglist));

%total parameters number
if stackmethod==3
  ncoef = nepoch; %multi-interferogram method
else
  ncoef = nifgs; %ifg by ifg method
end

dem = dem-dem(refy,refx);

for i=1:nifgs

  %coefficient for the linear atm correction
  if stackmethod==3
    %master/slave image epoch, using it to determine the cofficient position
    im = ifglist(i,1);
    is = ifglist(i,2);
    parm=atmparams(im);
    pars=atmparams(is);
    par = pars - parm;
  else
    par = atmparams(i);
  end

  %coefficients for the offsets
  joff = ncoef+i;

  %atm
  ifg_atm(:,:,i) = dem*par + atmparams(joff);
end
