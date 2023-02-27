function [Lon,Lat,interfero,R,dx,dy]=read_INTERFERO_tiff(file)
%function [Lon,Lat,interfero,R,dx,dy]=read_INTERFERO_tiff(rootname,track,volcano,master,slave)
% function to read the Interfero geotiff data
%track='079D_07694_131313';
%volcano='alayta';
%master='20161211';  
%slave='20170104';
 
% dirname=strcat(track,'/INTERF/');
% filename=strcat(volcano,'_',master,'_',slave,'.geo.unw.tif');
% file=strcat(rootname,dirname,filename);
 
[interfero,R] = geotiffread(file);
 
i=R.RasterSize(1);j=R.RasterSize(2);
Latmin = R.LatitudeLimits(1);
Latmax = R.LatitudeLimits(2);
Lonmin = R.LongitudeLimits(1);
Lonmax = R.LongitudeLimits(2);
dx = (Lonmax-Lonmin)./(j-1);
dy = (Latmax-Latmin)./(i-1);
 
Lon=Lonmin + dx.*[1:j] -dx;
Lat=Latmax - dy.*[1:i] +dy;
 
% cd(rootname)
% cd(dirname)
