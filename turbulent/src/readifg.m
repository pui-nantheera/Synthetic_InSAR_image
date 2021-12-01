%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  readifg.m: read interferogram              %
%%  date: 13/03/2008                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[ifg]=readifg(obsdir,ifg_nml,rows,cols) 

nifgs = length(ifg_nml);

% read interferograms 
for i=1:nifgs
  ifgname = char(strcat(obsdir,ifg_nml(i)));
  ifg_id = fopen(ifgname,'r','l');
  temp = fread(ifg_id, [cols,rows], 'real*4');
  fclose(ifg_id);
  ifg(:,:,i)=temp';
end
