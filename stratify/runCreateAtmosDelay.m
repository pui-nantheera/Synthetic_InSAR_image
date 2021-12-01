clear all
addpath('..\..\trieste_defo');
addpath('..\..\turbulent');

%% generate GACOS info
% -------------------------------------------------------------------------
rootname = 'G:\VolcanicUnrest\Atmosphere\stratify\volcano_2018\europe_etna\';

incidence=43.7835;
master = '20150805';
slave  = '20150817';

% reference point
xref  = 115.1;
yref  = -8.4;
dxref = 10;
dyref = dxref;

% parameters
wavelength = 0.055465;
m2rad = 4.*pi./wavelength;
rad2m = (wavelength./(4.*pi));
zen2los = 1./cos(incidence./180.*pi);

[~,~,atmo1] = read_GACOS([rootname,master]);
[~,~,atmo2] = read_GACOS([rootname,slave]);

atmo = (atmo2-atmo1).*zen2los.*m2rad;
atmo = imresize(atmo, [500 500]);

%% generate turbulent
% -------------------------------------------------------------------------
rows = 100;
cols = 100;
psizex = 1;
psizey = 1;
covmodel_type = 0;
maxvar = 7.5-1.5;       % mm
alpha = 1/12.3/10;     % km
N = 1;              % number of gaussian noise vectors
atm_pets = pcmc_atm(rows,cols,maxvar,alpha,covmodel_type,N,psizex,psizey);
atm_pets = imresize(atm_pets,[500 500]);

%% generate deformation
% -------------------------------------------------------------------------
Source_Type = 1;

% Source_Type = 1. Earthquakes
Quake.Strike = 25-15;              %strike in degrees 
Quake.Dip = 80-10;                 %dip in degrees
Quake.Rake = -90;               %rake in degrees
Quake.Slip = 1;                 %magnitude of slip vector in metres
Quake.Top_depth = 1+2.5 ;           %depth (measured vertically) to top of fault in kilometres
Quake.Bottom_depth = 10-5 ;       %depth (measured vertically) to bottom of fault in kilometres
Quake.Length = 10 -2;             %fault length in kilometres
% Source_Type = 2. Dykes
Dyke.Strike = 0-45;                %strike in degrees [0-180]
Dyke.Dip = 90-5;                  %dip in degrees (usually 90 or near 90)
Dyke.Opening = 1;               %magnitude of opening (perpendincular to plane) in metres
Dyke.Top_depth = 2-1.9 ;            %depth (measured vertically) to top of dyke in kilometres
Dyke.Bottom_depth = 8-2 ;         %depth (measured vertically) to bottom of dyke in kilometres
Dyke.Length = 10-5 ;              %dyke length in kilometres
% Source_Type = 3. Rectangular Sills
Sill.Strike = 0;                %strike (orientation of Length dimension) in degrees [no different]
Sill.Dip = 0;                   %Dip in degrees (usually zero or near zero)
Sill.Opening = 10;              %magnitude of opening (perpendincular to plane) in metres
Sill.Depth = 5;                 %depth (measured vertically) to top of dyke in kilometres
Sill.Width = 1;                 %depth (measured vertically) to bottom of dyke in kilometres
Sill.Length = 1;                %dyke length in kilometres
% Source_Type = 4. Magma Chamber - point pressure
Mogi.Depth  = 5;                %Depth of Mogi Source (5)
Mogi.Volume = 10*1e6;          	%Volume in m^3
% Source_Type = 5. Pressurized Penny-shaped Horizontal Crack (Fialko) - Sill
% Note, this is the slowest to calculate of the various sources
Penny.Depth  = 5;             	%Depth of crack in km^3
Penny.Pressure = 1*1e6;         %Pressure of crack in Pa
Penny.Radius  = 5;             	%Radius of crack in km^3
Heading = 192.04; 
x=-25000:100:25000-100;
y=-25000:100:25000-100;

[~, los_grid] = generateDeformation(Source_Type, x, y, Quake, Dyke, Sill, Mogi, Penny, Heading, incidence);
% scaling
los_grid = los_grid/0.028333*2*pi;

% combine all compoments
% -------------------------------------------------------------------------

weightS = 1;
weightT = 1;
sys_interfo = los_grid + weightS*atmo + weightT*atm_pets;
mask = imerode(atmo~=0,strel('disk',3));

% display result
% -------------------------------------------------------------------------

% %figure; set(gcf,'color','w');
% subplot(2,2,1); imagesc((wrapTo2Pi(los_grid)-pi).*mask); colormap jet; 
% axis image; title('deformation (D)');
% h=colorbar('vert');
% ylabel(h,'radians');
% subplot(2,2,2); imagesc((wrapTo2Pi(weightS*atmo)-pi).*mask); colormap jet; 
% axis image; title('stratified atmosphere (S)');
% h=colorbar('vert');
% ylabel(h,'radians');
% subplot(2,2,3); imagesc((wrapTo2Pi(weightT*atm_pets)-pi).*mask); colormap jet; 
% axis image; title('turbulent atmosphere (T)'); 
% h=colorbar('vert');
% ylabel(h,'radians');
% subplot(2,2,4); imagesc((wrapTo2Pi(sys_interfo)-pi).*mask); colormap jet; 
% axis image; title(['D + ',sprintf('%.2f',weightS),'S + ',sprintf('%.2f',weightT),'T']);
% h=colorbar('vert');
% ylabel(h,'radians');
% %%
% figure; set(gcf,'color','w');
% subplot(2,2,1); imagesc(los_grid.*mask); colormap jet; 
% axis image; title('deformation (D)');
% h=colorbar('vert');
% ylabel(h,'radians');
% subplot(2,2,2); imagesc((weightS*atmo).*mask); colormap jet; 
% axis image; title('stratified atmosphere (S)');
% h=colorbar('vert');
% ylabel(h,'radians');
% subplot(2,2,3); imagesc((weightT*atm_pets).*mask); colormap jet; 
% axis image; title('turbulent atmosphere (T)'); 
% h=colorbar('vert');
% ylabel(h,'radians');
% subplot(2,2,4); imagesc(sys_interfo.*mask); colormap jet; 
% axis image; title(['D + ',sprintf('%.2f',weightS),'S + ',sprintf('%.2f',weightT),'T']);
% h=colorbar('vert');
% ylabel(h,'radians');
% 
% %%
% 
% figure; set(gcf,'color','w');
subplot(2,4,2); imagesc((wrapTo2Pi(los_grid)-pi).*mask); colormap jet; 
axis image; title('wrapped deformation');
h=colorbar('vert');
ylabel(h,'radians');
subplot(2,4,4); imagesc((wrapTo2Pi(weightS*atmo)-pi).*mask); colormap jet; 
axis image; title('wrapped stratified atmosphere');
h=colorbar('vert');
ylabel(h,'radians');
subplot(2,4,6); imagesc((wrapTo2Pi(weightT*atm_pets)-pi).*mask); colormap jet; 
axis image; title('wrapped turbulent atmosphere'); 
h=colorbar('vert');
ylabel(h,'radians');
subplot(2,4,8); imagesc((wrapTo2Pi(sys_interfo)-pi).*mask); colormap jet; 
axis image; title(['wrapped D + S + T']);
h=colorbar('vert');
ylabel(h,'radians');

subplot(2,4,1); imagesc(los_grid.*mask); colormap jet; 
axis image; title('unwrapped deformation (D)');
h=colorbar('vert');
ylabel(h,'radians');
subplot(2,4,3); imagesc((weightS*atmo).*mask); colormap jet; 
axis image; title('unwrapped stratified atmosphere (S)');
h=colorbar('vert');
ylabel(h,'radians');
subplot(2,4,5); imagesc((weightT*atm_pets).*mask); colormap jet; 
axis image; title('unwrapped turbulent atmosphere (T)'); 
h=colorbar('vert');
ylabel(h,'radians');
subplot(2,4,7); imagesc(sys_interfo.*mask); colormap jet; 
axis image; title(['unwrapped D + S + T']);
h=colorbar('vert');
ylabel(h,'radians');