function [ifg_atm]=atmcorrect(ifg,ifglist,dem,lksx,lksy,refx,refy,stackmethod)

%% multilook interferogram
nifgs=length(ifglist);
for i=1:nifgs
  [ifg_lowres(:,:,i)]=looks(ifg(:,:,i),lksx,lksy);
end
[dem_lowres]=looks(dem,lksx,lksy);

[rows,cols,nifgs]=size(ifg);
rows_lowres=floor(rows/lksy);
cols_lowres=floor(cols/lksx);
refy_lowres=floor(refy/lksy);
refx_lowres=floor(refx/lksx);
n_lowres=rows_lowres*cols_lowres;

%%observations from matrix to vector
obsv=zeros(nifgs*n_lowres,1);
for i=1:nifgs
  i1 = (i-1)*n_lowres+1;
  i2 = i*n_lowres;
  obsv(i1:i2,1)=reshape(ifg_lowres(:,:,i)',n_lowres,1);
end

%%design matrix
[B]=atmdesign(ifglist,dem_lowres,refx_lowres,refy_lowres,stackmethod);

%delete NaN data points
B(isnan(obsv)==1,:)=[];
obsv(isnan(obsv)==1,:)=[];

%calculate atm parameters
if stackmethod==3
  atmparams=pinv(B)*obsv;
else
  atmparams=inv(B'*B)*(B'*obsv);
end
clear B;
clear obsv;

%%%design matrix for the full resolution interferogram
%forward calculation for the interferogram one by one considering the memory
[ifg_atm]=atmfwd(ifglist,dem,refx,refy,stackmethod,atmparams);
ifg_atm(isnan(ifg))=NaN;
