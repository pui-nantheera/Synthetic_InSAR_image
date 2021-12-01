function [Lon,Lat,interfero,dx,dy]=read_INTERFERO(track,master,slave)

% function to read the GACOS data
% track: 'ascending' or 'descending'
% master:  put the master acquisition
% slave:   put the slave acquisition

rootname=strcat('/users/fa17101/Documents/AGUNG/','INTERF','_',track,'/');
filename=strcat(master,'_',slave,'.diff_sm.unw.geo');
file=strcat(rootname,filename);

header = strcat(rootname,'EQA.dem_par');
%Data=importdata(header);
fid = fopen(header);
FC = textscan(fid, '%s%s%s%s');
fclose(fid);

xi=str2num(cell2mat(FC{2}(11)));
yi=str2num(cell2mat(FC{2}(10)));
dx=str2num(cell2mat(FC{2}(13)));
dy=-dx;
j=str2num(cell2mat(FC{2}(8)));
i=str2num(cell2mat(FC{2}(9)));

% if strcmp(master,'20170921')==1
%     i=i-1;
% end

interfero = multibandread(file, [i j 1], ...
                     'float', 0, 'bsq', 'ieee-be', ...
                     {'Row', 'Range', [1 i]}, ...
                    {'Column', 'Range', [1 j]} );

Lon=xi + dx.*[1:j] - dx./2;
Lat=yi + dy.*[1:i] - dy./2;