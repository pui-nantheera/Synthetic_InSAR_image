clear all
close all
clc

% interfero
track='descending';
incidence=43.7835;
master='20170828';
slave ='20170909';

% reference point
xref=115.1;
yref=-8.4;
dxref=10;
dyref=dxref;

% parameters
wavelength = 0.055465;
m2rad=4.*pi./wavelength;
rad2m=(wavelength./(4.*pi));
zen2los=1./cos(incidence./180.*pi);

[Lon1,Lat1,atmo1,dx1,dy1]=read_GACOS(track,master);
[Lon1,Lat1,atmo2,dx1,dy1]=read_GACOS(track,slave);
atmo=((atmo2-atmo1).*(zen2los).*(m2rad));

[Lon2,Lat2,inter,dx2,dy2]=read_INTERFERO(track,master,slave);


%% resize interfero same dx,dy than atmo
[i,j]=size(inter);
i2=(i./round(dx1./dx2));
j2=(j./round(dx1./dx2));
INTER = imresize(inter,[i2 j2]);
Lon=Lon2(1) + dx1.*[1:j2] - dx1./2;
Lat=Lat2(1) + dy1.*[1:i2] - dy1./2;

%% crop atmo to have same size at interfero
Loni=Lon(1);Lonf=Lon(end);
Lati=Lat(1);Latf=Lat(end);
imin=find(abs(Lat1-Lati)==min(abs(Lat1-Lati)));
imax=find(abs(Lat1-Latf)==min(abs(Lat1-Latf)));
jmin=find(abs(Lon1-Loni)==min(abs(Lon1-Loni)));
jmax=find(abs(Lon1-Lonf)==min(abs(Lon1-Lonf)));

ATMO=atmo(imin:imax,jmin:jmax);
LAT=Lat1(imin:imax);
LON=Lon1(jmin:jmax);

%% mask zero values
MASK=zeros(i2,j2);
MASK(INTER~=0)=1;

%% difference from ref
iref=find(abs(LAT-yref)==min(abs(LAT-yref)));
jref=find(abs(LON-xref)==min(abs(LON-xref)));
sub1=INTER(iref-dxref:iref+dxref,jref-dyref:jref+dyref);
DS1=sub1(sub1~=0);
INTER_ref=mean(DS1);
sub2=ATMO(iref-dxref:iref+dxref,jref-dyref:jref+dyref);
DS2=sub2(sub2~=0);
ATMO_ref=mean(DS2);
INTER=(INTER-INTER_ref).*MASK;
ATMO=(ATMO-ATMO_ref).*MASK;

%% figure
figure(1)
subplot(2,2,1)
imagesc(LON,LAT,INTER);
title(strcat('Interferogram-',master,'-',slave));
colormap jet
colorbar
axis xy
caxis([-15 15])

subplot(2,2,2)
imagesc(LON,LAT,ATMO);
title(strcat('GACOS model-',master,'-',slave));
colormap jet
colorbar
axis xy
caxis([-15 15])

subplot(2,2,3)
imagesc(LON,LAT,INTER-ATMO);
title(strcat('Interferogram-corrected-',master,'-',slave));
colormap jet
colorbar
axis xy
caxis([-15 15])