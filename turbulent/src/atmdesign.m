function [B]=atmdesign(ifglist,dem,refx,refy,stackmethod)

nifgs=length(ifglist);
[rows,cols]=size(dem);
nepoch = max(max(ifglist));
n = rows*cols;

%total parameters number
if stackmethod==3
  ncoef = nepoch; %multi-interferogram method
else
  ncoef = nifgs; %ifg by ifg method
end
npar = ncoef+nifgs; %total parameters number

dem = dem-dem(refy,refx);
demv = reshape(dem',n,1);

B = zeros(n*nifgs,npar);
for i=1:nifgs
  %master/slave image epoch, using it to determine the cofficient position
  im = ifglist(i,1);
  is = ifglist(i,2);
  %start/end line number for each interferogram in the coefficient matrix
  ib1 = (i-1)*n+1;
  ib2 = i*n;

  %coefficient for the linear atm correction
  if stackmethod==3
    B(ib1:ib2,im)=-demv;
    B(ib1:ib2,is)=demv;
  else
    B(ib1:ib2,i)=demv;
  end

  %coefficients for the offsets
  joff = ncoef+i;
  B(ib1:ib2,joff)=ones(n,1);
end
