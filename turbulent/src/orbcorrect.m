function [ifg_orb,tiltpar]=orbcorrect(ifg,ifglist,order,lksx,lksy,refx,refy,psize,stackmethod)

%% multilook interferogram
[rows,cols,nifgs]=size(ifg);
for i=1:nifgs
  [ifg_lowres(:,:,i)]=looks(ifg(:,:,i),lksx,lksy);
end
rows_lowres=floor(rows/lksy);
cols_lowres=floor(cols/lksx);
n_lowres=rows_lowres*cols_lowres;
refy_lowres=floor(refy/lksy);
refx_lowres=floor(refx/lksx);
psizey_lowres=psize*lksy;
psizex_lowres=psize*lksx;

%%observations from matrix to vector
obsv=zeros(nifgs*n_lowres,1);
for i=1:nifgs
  i1 = (i-1)*n_lowres+1;
  i2 = i*n_lowres;
  obsv(i1:i2,1)=reshape(ifg_lowres(:,:,i)',n_lowres,1);
end
clear ifg_lowres;

%%design matrix
[B]=orbdesign(ifglist,order,rows_lowres,cols_lowres,refx_lowres,refy_lowres,psizex_lowres,psizey_lowres,stackmethod);

%delete NaN data points
B(isnan(obsv)==1,:)=[]; 
obsv(isnan(obsv)==1,:)=[]; 

%calculate orbit parameters
if stackmethod==3
  orbparams=pinv(B)*obsv;
else
  invNbb = inv(B'*B);
  orbparams=invNbb*(B'*obsv);
end
clear B;
clear obsv;

%%%orbitdesign for the full resolution interferogram
%forward calculation for the interferogram one by one considering the memory
[ifg_orb] = orbfwd(ifglist,order,rows,cols,refx,refy,psize,psize,stackmethod,orbparams);
ifg_orb(isnan(ifg))=NaN;

% reshape the parameter and output
%calculate orbit parameters, without the constant offset
ncoef = (order+1)*(order+2)/2-1;
ncoefline = floor((length(orbparams)-nifgs)/ncoef);
tiltpar = reshape(orbparams(1:ncoefline*ncoef),ncoef,ncoefline);
%stdtiltpar = std(tiltpar');
