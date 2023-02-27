function [Lon,Lat,atmo,dx,dy]=read_GACOS_tiff(file, header)
 
% function to read the Interfero geotiff data
%track='079D_07694_131313';
%volcano='alayta';
%master='20170104';
 
%dirname=strcat(track,'/GACOS/');
%file = strcat(rootname,dirname,master,'.ztd');
%header = strcat(file,'.ztd.rsc');
fid = fopen(header);
FC = textscan(fid, '%s%s%s%s');
fclose(fid);
 
xi=str2num(cell2mat(FC{2}(7)));
yi=str2num(cell2mat(FC{2}(8)));
dx=str2num(cell2mat(FC{2}(9)));
dy=str2num(cell2mat(FC{2}(10)));
j=str2num(cell2mat(FC{2}(1)));
i=str2num(cell2mat(FC{2}(2)));
 
atmo = multibandread(file, [i j 1], ...
                     'float', 0, 'bsq', 'ieee-le', ...
                     {'Row', 'Range', [1 i]}, ...
                    {'Column', 'Range', [1 j]} );
                
Lon=xi + dx.*[1:j] - dx./2;
Lat=yi + dy.*[1:i] - dy./2;
