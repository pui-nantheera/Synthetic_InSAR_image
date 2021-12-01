%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  tsinv.m: time series inversion pixel by pixel    %
%%  author: Hua Wang                                 %
%%  date: 27/03/2008                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ifg_ts]=tsinv(ifg,ifglist,epochlist,vcm,pthr)

[rows,cols,nifgs]=size(ifg);

nepoch=max(max(ifglist));
npar = nepoch-1;

%% fill in the design matrix
B0 = zeros(nifgs,npar);
for i=1:nifgs
  im = ifglist(i,1);
  is = ifglist(i,2);
  if is<im
    tmp=im;
    im=is;    
    is=tmp;
  end
  B0(i,im:is-1)=epochlist(im+1:is,3);
end

%% calculate the time series
ifg_ts=NaN(rows,cols,npar);
ifg = permute(ifg,[3,1,2]);

for i=1:rows
  for j=1:cols
    %initialize
    B=B0;
    ifgv = ifg(:,i,j);
    %delete NaN pixels
    B(isnan(ifgv),:)=[];
    m=size(B,1);

    %calculate slip rate by iterative least-squares
    if (m>=pthr)
      vcm_tmp = vcm;
      vcm_tmp(isnan(ifgv),:)=[];
      vcm_tmp(:,isnan(ifgv))=[];
      ifgv(isnan(ifgv)==1)=[];
      %add the laplacian observation
      P = pinv(vcm_tmp);  %weight matrix
      Nbb=B'*P*B;
      W = B'*P*ifgv;
      x = pinv(Nbb)*W;  %slip rate on pixel for each epoch
      ifg_ts(i,j,:)=x;
    else
      ifg_ts(i,j,:)=NaN;
    end
  end
end

for i=1:npar
  ifg_ts(:,:,i)=ifg_ts(:,:,i)*epochlist(i+1,3);
end
