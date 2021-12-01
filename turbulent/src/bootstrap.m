%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bootstrap.m: program for bootstrap method                 %
%%  date: 17/04/2008                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = bootstrap()

clear all

%%%% input interferograms directory %%%%
ifg_list = 'ifgall.list';
ifgsub_list = 'ifg.list';
root   = '/nfs/see-archive-01_a1/earhw/xsh/analysis/';
obsdir = strcat(root,'obs/');

%number of interferograms to be removed
nrm = 2;
%datasets for bootstrap
nsets=100;

%%% sort the master/slave namelist here
ifgfilelist = strcat(obsdir,ifg_list);
ifgsubfilelist = strcat(obsdir,ifgsub_list);
[ifg_nml] = textread(ifgfilelist,'%s');

nifgs = length(ifg_nml);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('---------bootstrap method---------------');
for i=1:nrm
  irm(:,i) = ceil((nifgs-i+1)*rand(nsets,1));
end

for iset=1:nsets
  display(strcat('running...',num2str(iset)));
  display(strcat('remove interferogram:',' ( ',char(ifg_nml(irm(iset,:))),')'));
  ifgsub_nml = ifg_nml;
  for i=1:nrm
    ifgsub_nml(irm(iset,i)) = [];
  end
  fid = fopen(ifgsubfilelist,'w');
  for i=1:nifgs-nrm
    fprintf(fid,'%30s\n',char(ifgsub_nml(i)));
  end
  fclose(fid);
  sliprateset(iset) = main();
  display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
end

sliprateset(isnan(sliprateset))=[]; %delete non-convergent solutions
finalsliprate = mean(sliprateset);
errorsliprate = std(sliprateset);
display(num2str(sliprateset));
display(strcat('final slip rate:',num2str(finalsliprate),'+-',num2str(errorsliprate)));
display(strcat('maximum slip rate:',num2str(max(sliprateset))));
display(strcat('minimum slip rate:',num2str(min(sliprateset))));
display(strcat('final convergent number:',num2str(length(sliprateset))));
