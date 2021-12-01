%
%getparm: get a parameter from a two-column matrix
%
function [parval] = getpar(parname,parmat,fmt)

npar = size(parmat,1);
for i=1:npar
  strtmp = parmat(i,1);
  if strcmp(parname,strtmp)==1
    parval = parmat(i,2);
    break;
  end
end

if i>npar
  display(['Can not find parameter: ' parname]);
  exit
end

parval=char(parval);
if strcmp(fmt,'n')==1
  parval=str2num(char(parval));
end
