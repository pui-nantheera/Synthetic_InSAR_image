%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  jackknife.m: program for jackknife method                 %
%%  date: 17/04/2008                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = jackknife()

clear all

%%% input interferograms directory %%%%
ifg_list = 'ifgall.list';
ifgsub_list = 'ifg.list';
root   = '/nfs/see-archive-01_a1/earhw/xsh/analysis/';
obsdir = strcat(root,'obs/');

%%% sort the master/slave namelist here
ifgfilelist = strcat(obsdir,ifg_list);
ifgsubfilelist = strcat(obsdir,ifgsub_list);
[ifg_nml] = textread(ifgfilelist,'%s');

%%% jack slip rate inversion
nifgs = length(ifg_nml);
for i=1:nifgs
  display(strcat('remove interferogram:',num2str(i),' ( ',char(ifg_nml(i)),')'));
  ifgsub_nml = ifg_nml;
  ifgsub_nml(i) = [];
  fid = fopen(ifgsubfilelist,'w');
  for j=1:nifgs-1
    fprintf(fid,'%30s\n',char(ifgsub_nml(j)));
  end
  fclose(fid);
  sliprate(i) = main();
  display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
end

jack_sliprate = sliprate;
jack_sliprate(isnan(sliprate))=[]; %delete non-convergent solutions
jack_finalsliprate = mean(jack_sliprate);
jack_errorsliprate = std(jack_sliprate);
display(num2str(sliprate));
display(strcat('final slip rate (jackknife method):',num2str(jack_finalsliprate),'+-',num2str(jack_errorsliprate)));
display(strcat('maximum slip rate (jackknife method):',num2str(max(jack_sliprate))));
display(strcat('minimum slip rate (jackknife method):',num2str(min(jack_sliprate))));
display(strcat('final convergent number (jackknife method):',num2str(length(jack_sliprate))));
