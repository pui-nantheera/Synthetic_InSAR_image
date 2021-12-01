%
%readparfile: get parameter matrix from a configure file
%
function [parmat] = readparfile(parfile)

parmat=[];
fid = fopen(parfile,'r');
while (feof(fid)==0)
  strline = fgetl(fid);
  if length(strline)>1
    if (strline(1)~='%' && strline(1)~='#')
      [parname,parval]=strread(strline,'%s %s');
      parmat=[parmat;parname,parval];
    end
  end
end
fclose(fid);
