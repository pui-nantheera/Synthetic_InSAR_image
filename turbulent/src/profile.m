%%%% Author: Hua Wang, 04 Mar 2008

function [profile,sumprof] = profile(ifg,pfault,swath,nwin,psize,stdifg,alpha)

%pfault: the coordinates of the profile perpendicular to the fault
[rows,cols]=size(ifg);
n=rows*cols;

if nargin<6
  stdifg = ones(rows,cols);
end

%get the transform matrix
%[ cos(alpha)  sin(alpha) ]
%[-sin(alpha)  cos(alpha) ]
sx = pfault(2,1)-pfault(1,1);  %x coordinate
sy = pfault(2,2)-pfault(1,2);  %y coordinate
l = sqrt(sx*sx+sy*sy);
coef = [sx/l, sy/l; -sy/l, sx/l];

%coordinate transform
[yy,xx]=meshgrid(1:rows,1:cols);
xxv=reshape(xx,n,1);
yyv=reshape(yy,n,1);
grid=[xxv-pfault(1,1), yyv-pfault(1,2)];
grid=grid';
outgrid = coef*grid;
outgrid = (outgrid').*psize; %convert to km
clear xxv;
clear yyv;
clear xx;
clear yy;

%dataset including outgrid,ifg,stdifg
ifgv=reshape(ifg',n,1);
stdv=reshape(stdifg',n,1);

%profile(:,1): distance along the profile
%profile(:,2): distance normal to the profile
%profile(:,3): ifg
%profile(:,4): std of ifg
%remove pixels outside of swath and NAN
profile=[outgrid,ifgv,stdv];
profile(isnan(profile(:,3)),:)=[];
profile=profile(find(profile(:,1)>0 & profile(:,1)<l*psize & abs(profile(:,2))<swath),:);

clear ifg;
clear ifgstd;
clear outgrid;

k=1;
smin = floor(min(profile(:,1)));
smax = floor(max(profile(:,1)));
for ibin=smin:nwin:smax-nwin
  bin=profile(find(profile(:,1)>=ibin & profile(:,1)<(ibin+nwin)),:); 
  npix = size(bin,1);
  if npix>1
    %variance-covariance matrix
    if nargin>5
      xdistdiff=repmat(bin(:,1),1,npix)-repmat(bin(:,1)',npix,1);
      ydistdiff=repmat(bin(:,2),1,npix)-repmat(bin(:,2)',npix,1);
      dist=sqrt(xdistdiff.^2+ydistdiff.^2);
      vcm=exp(-dist.*alpha);
      vcm_tmp = bin(:,4)*bin(:,4)';
      vcm = vcm_tmp.*vcm;
    else
      vcm=diag(bin(:,4));
    end
    %ls inversion
    B = ones(npix,1);
    [x,stdx]=lscov(B,bin(:,3),vcm);
    distprof=mean(bin(:,1));
    sumprof(k,:)=[distprof,x,stdx];
    k=k+1;
  end
end
