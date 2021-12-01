%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% make_vcm: calculate the variance-covariance matrix  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vcm_t,vcm_s,maxvar,alpha] = make_vcm(ifg_flat,ifglist,psize,lksx,lksy)

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

%%% set up template for variance-covariance matrix. c=0.5 for common master
%%% or slave; c=-0.5 if master of one matches slave of another
vcm_pat=diag(ones(nifgs,1));

disp('calculating vcm in time domain ...');
for i = 1:nifgs-1
  %master/slave image name
  mas1 = ifglist(i,1);
  slv1 = ifglist(i,2);
  for j = i+1:nifgs
    %master/slave image name
    mas2 = ifglist(j,1);
    slv2 = ifglist(j,2);
    if (mas1==mas2 || slv1==slv2)
      vcm_pat(i,j) = 0.5;
    end
    if (mas1==slv2 || slv1==mas2)
      vcm_pat(i,j) = -0.5;
    end
    vcm_pat(j,i)=vcm_pat(i,j);
  end
end

%%% make covariance matrix in time domain
var = sqrt(maxvar);
vcm_t = var*var';
vcm_t = vcm_t.*vcm_pat;

%%% make the covariance matrix in the spatial domain
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
