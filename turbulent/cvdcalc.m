function [maxvar,alpha] = cvdcalc(flatgrid,ncols,nlines,pixsize,method,saved_acg)
% function [cvdav,maxvar,alpha,acg,expcosmod,bessmod] = cvdcalc(flatgrid,ncols,nlines,pixsize,method,saved_acg)
%function [cvdav,maxvar,alpha,acg,expcosmod,bessmod,flatgrid] = cvdcalc(grid,ncols,nlines,pixsize,nullval,method,byteorder,saved_acg)
% function [cvdav,maxvar,alpha,acg,expcosmod,bessmod,flatgrid] = cvdcalc(infile,ncols,nlines,pixsize,nullval,method,byteorder)
%
% or...
%
% cvdcalc(...,saved_acg)   if have a autocorr_grid file saved as a .mat file [variable must be named autocorr_grid] 
%
% Calculate average covariance vs distance (autocorrelation) function
% and best-fitting exponential to it.
%
% Does calculation spatially (slow) using xcorr2, or spectrally using Wiener-Khinchine Theorem
%
% INPUTS:
%
% infile = filename of binary grid (assume real*4, big-endian (MSBfirst))
% ncols  = number of columns
% nlines = number of lines
% pixsize = pixel size in kilometres
% nullval = null value of grid
% method = 1 (do calculation spectrally) or 2 (do calc spatially)  
% byteorder = 'b'/'l' for big endian (MSBFirst) / little endian (LSBFirst)
% saved_acg = 'save_acg'.mat file containing variable autocorr_grid from a previous run of cvdcalc
%
% OUTPUTS:
%
% cvdav = 2 column file. col1 = distance, col2 = covariance
% maxvar = maximum covariance
% alpha = bestfit exponential decay model (covariance = maxvar(-alpha*dist))
% acg = autocorrelation grid
% expcosmod = parameters b and c from maxvar*exp(-br)*cos(cr)
% bessmod = parameters r and w from maxvar*exp(-x/r)*J0(2*pi*x/w)
%
  % requires flatten.m, expcos.m, pendiffexpcos.m, ebessel.m, pendiffebessel.m
%
% tjw 17-jul-02
%

% modifications
% added spectral method, 11 feb 2003

%clf

%%display grid and flattened grid
%sort out colortable
% colormap('default')
% map = colormap;
% map(1,:)=[0 0 0];
% colormap(map);
%%%display grid
% subplot(2,2,1)
% dispgrid = flatgrid;
% dispgrid(flatgrid==0) = min(min(flatgrid))-0.00001;
% image(dispgrid,'CDataMapping','scaled')
% colorbar('vert')
% title('Input interferogram')
%%%display flattened grid
% subplot(2,2,2)
% dispgrid = flatgrid;
% dispgrid(flatgrid==0) = min(min(flatgrid))-0.00001;
% image(dispgrid,'CDataMapping','scaled')
% colorbar('vert')
% title('Flattened interferogram')
 
%% calculate autocorrelation grid
%disp(['calculating autocorrelation grid...'])
if nargin<8
  if method==2
      autocorr_grid = xcorr2(flatgrid);
      autocorr_grid = autocorr_grid./nnz(flatgrid);
  elseif method==1
      fft_flatgrid=fft2(flatgrid);
      pspec = real(fft_flatgrid).^2 + imag(fft_flatgrid).^2;
      autocorr_grid = ifft2(pspec);
      autocorr_grid = fftshift(real(autocorr_grid))/nnz(flatgrid);         
  else
      disp(['Method flag in input must = 1 (spectrally) or 2 (spatially)'])
      exit
  end
else 
    load(saved_acg)
end

%disp(['...finished'])
if method==2
  xcent = ncols;
  ycent = nlines;
  ncols =         ncols*2 -1 ; %ncols/nlines now = of autocorr_grid
  nlines =        nlines*2-1 ;
else
  xcent = floor(ncols/2)+1 ;   % spectrally, the grid is not doubled in size
  ycent = floor(nlines/2)+1 ;
end

%%% display autocorr_grid
%subplot(2,2,3)
%imagesc (autocorr_grid)
%colorbar('vert')
%title('Autocorrelation grid')

%% make grid of distances from centre 
[xx,yy] = meshgrid(1:ncols,1:nlines);
r = sqrt((xx-xcent).^2+(yy-ycent).^2)*pixsize;
r = reshape(r,ncols*nlines,1);
acg = reshape(autocorr_grid,ncols*nlines,1);

if method==2
  r = r(1:ceil(length(r)/2)); %only need 1st half of image (symmetry)
  acg= acg(1:length(r));
else
  r = r(1:(ceil(length(r)/2)+nlines)); %only need 1st half of image (symmetry)
  acg= acg(1:length(r));
end

%% fit exponential decay cov = maxcov*exp(-alpha*r) 
rorig = r; %keep copies of original columns
acgorig = acg;
%r(acg<=0)=[]; %blank out entries with covariance <=0
%acg(acg<=0)=[];
%acg(r>(xcent*pixsize/2))=[]; %blank out distant points
%r(r>(xcent*pixsize/2))=[];
maxacg=max(acg);
%r(acg<(maxacg/exp(1)))=[]; %blank out points where acg < maxacg/e
%acg(acg<(maxacg/exp(1)))=[];
%params = r\(log(acg)-log(maxacg)); % solve for alpha using least squares
%alpha = -params(1)
maxr = ceil(max(rorig));

%%calculate average cov vs dist profile
w = pixsize*2 ;  %bin width
if (xcent<ycent)
    maxdist = xcent*pixsize;
else
    maxdist = ycent*pixsize;
end
r=rorig;
r(rorig>=maxdist)=[];
length(r);
acg=acgorig;
acg(rorig>=maxdist)=[];
length(acg);
rbin=ceil(r./w); % classify values of r according to bin number
cvdav = zeros(max(rbin),2);
maxbin = max(rbin)-1;
for bin = 0:maxbin,
    cvdav(bin+1,2) = mean(acg(find(rbin==bin)));
    cvdav(bin+1,1) = bin*w;
end

%% calculate best fit function maxvar*exp(-alpha*r)
alpha = fminsearch('pendiffexp',[2/(maxbin*w)],[],cvdav);

%% calculate best fit function maxvar*exp(-br)*cos(cr)
%% expcosmod = fminsearch('pendiffexpcos',[alpha alpha]',[],cvdav)

%% calculate best fit function maxvar*exp(-x/r)*J0(2*pi*x/w)
%% bessmod = fminsearch('pendiffebessel',[exp(1)/alpha exp(1)/alpha]',[],cvdav)


%% plot cvd data, average and exponential fit
%subplot(2,2,4)
%plot(rorig,acgorig,'.')
%hold on
%plot(cvdav(:,1),cvdav(:,2),'g')
%plot(0:maxr,maxacg*exp(-alpha*(0:maxr)),'r')
%plot(cvdav(:,1),expcos(cvdav(:,1),maxacg,expcosmod(1),expcosmod(2)),'y')
%plot(cvdav(:,1),ebessel(cvdav(:,1),maxacg,bessmod(1),bessmod(2)),'m')

%%tidy up output
% acg = autocorr_grid;
maxvar=maxacg;
%save allparams.mat
