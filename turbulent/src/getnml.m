function [ifg_nml,BaseT,ifglist,epochlist] = getnml(ifgfilelist)

%read the interferogram file list
[ifg_nml] = textread(ifgfilelist,'%s');

%% get master and slave images from the interferogram namelist
nifgs=length(ifg_nml);
mas=zeros(nifgs,1);
slv=zeros(nifgs,1);
for i=1:nifgs;
  strpair=char(ifg_nml(i));
  mas(i) = str2double(strpair(5:10));
  slv(i) = str2double(strpair(12:17));
  %solve for the 2000 problem
  if mas(i)>500000
    mas(i) = mas(i)+19000000;
  else
    mas(i) = mas(i)+20000000;
  end
  if slv(i)>500000
    slv(i) = slv(i)+19000000;
  else
    slv(i) = slv(i)+20000000;
  end
end

%%% sort the sar images, including all the master & slave images
sarlist = sort([mas;slv]);

%%% incorporate the same images
%%% epochlist(:,1): epochlist in numeric format, e.g., 19991212
%%% epochlist(:,2): times used for this epoch, e.g., 2
%%% epochlist(:,3): epochlist in year format w.r.t the first epoch
nimages = length(sarlist);
epochlist(1,1) = sarlist(1);
epochlist(1,2) = 1;
j=2;
for i=2:nimages
  % new epoch
  if sarlist(i)~=sarlist(i-1)
    epochlist(j,1)=sarlist(i);
    epochlist(j,2)=1;
    j=j+1;
  else
    epochlist(j-1,2)=epochlist(j-1,2)+1;
  end
end
nimages=length(epochlist);

%convert from yyyymmdd2year
for i=1:nimages
  epochlist(i,3)=yyyymmdd2year(epochlist(i,1));
end
startdate = epochlist(1,3);
for i=1:nimages
  epochlist(i,3)=epochlist(i,3)-startdate;
end

%%% get the interferogram list
for i=1:nifgs
  for j=1:nimages
    if mas(i)==epochlist(j,1)
      ifglist(i,1)=j;
    end
    if slv(i)==epochlist(j,1)
      ifglist(i,2)=j;
    end
  end
end

%temporal baseline
BaseT = epochlist(ifglist(:,2),3) - epochlist(ifglist(:,1),3);
