%makes model interferograms based on sliprate estimate. 

function [fmodel]=make_model(fmodelfile,rows,cols)

%read in forward model
%NOTE: pay attention to the file format (real*4) and unit (mm/yr)
fid = fopen (fmodelfile,'r','l');
fmodel =fread(fid,[cols,rows],'real*4');
fclose(fid);

%comment by Hua Wang
%output is in mm/yr already using hua's cehm program
%fmodel=fmodel'.*10; %convert to mm/yr.
fmodel=fmodel';
