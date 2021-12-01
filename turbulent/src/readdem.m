%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  readdem.m: read dem data                   %
%%  date: 29/04/2008                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[dem]=readdem(demfile,rows,cols) 
dem_id = fopen(demfile,'r','l');
temp = fread(dem_id, [cols,rows], 'int16');
fclose(dem_id);
dem=temp';
