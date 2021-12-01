%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% make_vcm_correct: calculate the variance-covariance matrix  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vcm_t,vcm_s] = make_vcm_correct(ifg_flat,ifglist,epochlist,psize,lksx,lksy)

%%vcm_t: variance-covariance matrix in time domain
%%vcm_s: variance-covariance matrix in space domain

[rows,cols,nifgs]=size(ifg_flat);

disp('inverting atmospheric delay parameters ...');
ifg_flat(isnan(ifg_flat))=0;
maxvar=zeros(nifgs,1);
alpha=zeros(nifgs,1);
for i=1:nifgs
  [maxvar(i),alpha(i)] = cvdcalc(ifg_flat(:,:,i),cols,rows,psize,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% make the covariance matrix in time domain  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('calculating vcm in time domain ...');

% calculate maximum variance for each epoch using pinv method
nepoch = max(max(ifglist));

B = zeros(nifgs,nepoch);
for i=1:nifgs
  %master/slave image epoch, using it to determine the cofficient position
  im = ifglist(i,1);
  is = ifglist(i,2);
  %fill in the design matrix
  B(i,im)=1;
  B(i,is)=1;
end

varobs = maxvar-1;

% using the svd method 
%varepoch = pinv(B)*varobs;  

varepoch = lsqnonneg(B,varobs);
varepoch = varepoch+0.5;

% using the constrained linear least-squares method
% lb=ones(nepoch,1);  %supposing the variance is larger than 1 mm
% varepoch = lsqlin(B,varobs,[],[],[],[],lb,[]);

% generate the covariance matrix
vcm_t=diag(maxvar);
for i = 1:nifgs-1
  %master/slave image name
  mas1 = ifglist(i,1);
  slv1 = ifglist(i,2);
  for j = i+1:nifgs
    %master/slave image name
    mas2 = ifglist(j,1);
    slv2 = ifglist(j,2);
    if mas1==mas2
      vcm_t(i,j) = varepoch(mas1);
    elseif slv1==slv2
      vcm_t(i,j) = varepoch(slv1);
    elseif mas1==slv2
      vcm_t(i,j) = -varepoch(mas1);
    elseif mas2==slv1
      vcm_t(i,j) = -varepoch(mas2);
    end
    vcm_t(j,i)=vcm_t(i,j);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% make the covariance matrix in space domain %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('calculating vcm in space domain ...');

ma = mean(alpha);
rows_lowres=floor(rows/lksy);
cols_lowres=floor(cols/lksx);
psizex=lksx*psize;
psizey=lksy*psize;

[ydist,xdist]=meshgrid(1:rows_lowres,1:cols_lowres);
xdistv=reshape(xdist,rows_lowres*cols_lowres,1);
ydistv=reshape(ydist,rows_lowres*cols_lowres,1);
xdistv=xdistv.*psizex;
ydistv=ydistv.*psizey;
xdistdiff=repmat(xdistv,1,length(xdistv))-repmat(xdistv',length(xdistv),1);
ydistdiff=repmat(ydistv,1,length(xdistv))-repmat(ydistv',length(xdistv),1);
dist=sqrt(xdistdiff.^2+ydistdiff.^2);
vcm_s=exp(-dist.*ma);
