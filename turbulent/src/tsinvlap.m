%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  tsinv.m: time series inversion pixel by pixel    %
%%           with Laplacian smoothing                %
%%  author: Hua Wang                                 %
%%  date: 27/03/2008                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ifg_ts,ifg_fit,ifg_res,ifg_var,ts_rough]=tsinvlap(ifg,ifglist,epochlist,vcm,pthr,r2)

[rows,cols,nifgs]=size(ifg);

nepoch=max(max(ifglist));
npar = nepoch-1;
nlap = npar-2;
nobs=nifgs+nlap;

%% fill in the design matrix
B0 = zeros(nifgs,npar);
for i=1:nifgs
  im = ifglist(i,1);
  is = ifglist(i,2);
  if is<im
    min=im;
    im=is;    
    is=min;
  end  
  B0(i,im:is-1)=epochlist(im+1:is,3);
end

%% Laplacian smooth
% coefficient matrix
BLap = zeros(nlap,npar);
for i=1:nlap
  BLap(i,i) = 0.5;
  BLap(i,i+2) = 0.5;
  BLap(i,i+1) = -1;
end
% observation vector
vLap = zeros(nlap,1);

% total coefficient matrix
B0 = [B0;BLap];

%% new covariance matrix, adding the laplacian equations
vcm_new = diag(ones(nobs,1).*r2);
vcm_new(1:nifgs,1:nifgs)=vcm;

%% calculate the time series
ifg_ts=NaN(rows,cols,npar);
ifg_fit=NaN(rows,cols,nifgs);
ifg_var=NaN(rows,cols);
ts_rough=NaN(rows,cols);

ifg = permute(ifg,[3,1,2]);

for i=1:rows
  for j=1:cols
    %initialize
    ifgv = ifg(:,i,j);
    mask=ones(nifgs,1);
    mask(isnan(ifgv))=0;
    m=sum(mask);
    %calculate slip rate by iterative least-squares
    if (m>=pthr)
      B=B0;
      vcm_tmp = vcm_new;
      B(isnan(ifgv),:)=[];
      vcm_tmp(isnan(ifgv),:)=[];
      vcm_tmp(:,isnan(ifgv))=[];
      ifgv(isnan(ifgv)==1)=[];
      %add the laplacian observation
      ifgv = [ifgv;vLap];

      [x] = lscov(B,ifgv,vcm_tmp);
      fit = B0*x;
      ifg_fit(i,j,:)=fit(1:nifgs);
      ifg_ts(i,j,:)=x;

      %residuals of the fitting
      v = B*x-ifgv;
      ifg_var(i,j)=sqrt(v(1:m)'*v(1:m)/m);
      ts_rough(i,j)=sqrt(fit(m+1:)'*fit(m+1:)/nlap);
    end
  end
end

%calculate the final residuals
ifg_res=ifg-ifg_fit;

%accumulate time series
ifg_ts(:,:,1)=ifg_ts(:,:,1)*epochlist(2,3);
for i=2:npar
  ifg_ts(:,:,i)=ifg_ts(:,:,i-1)+ifg_ts(:,:,i)*epochlist(i+1,3);
end
