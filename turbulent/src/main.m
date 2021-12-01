%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  main.m: main program for multi-interferogram method       %
%%  date: 13/03/2008                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [finalsliprate] = main(parfile)

 function [finalsliprate] = main()
 parfile='10km_K0.5delta0.conf';

%read parameter file
[parmat] = readparfile(parfile);

%%%% interferogram size %%%%%%%%
rows  = getpar('rows:',parmat,'n');
cols  = getpar('cols:',parmat,'n');
psize = getpar('psize:',parmat,'n');

%%% Option to solve for two components of slip 
%%% horzontal and vertical(rake='variable', two input forward models)
%%% or a single component (rake='fixed', single input forward model)
rake = getpar('rake:',parmat,'s');

%%%% slip rate for the initial forward model in mm/yr
inisliprate_horz = getpar('inisliprate_horz:',parmat,'n');
%%%% slip rate for the synthetic data mm/yr, used to find the folder of the synthetic data
simsliprate_horz = getpar('simsliprate_horz:',parmat,'n');

%%% number of sets of data (multiple sets required for Monte Carlo error estimates %%%
%realdata: enter 0 for synthetic/montecarlo data, 1 for a single set of real data.
realdata = getpar('realdata:',parmat,'n');
if realdata==0
  nsets = getpar('nsets:',parmat,'n');
else
  nsets = 1;
end

%%% stacking parameters
%stackmethod 1: classic stacking; 2: pixel by pixel stacking; 3: multi-interferogram method
stackmethod = getpar('stackmethod:',parmat,'n');
%nsig: for temporal low-pass filter (sigma < nsig*sigma0)
nsig = getpar('nsig:',parmat,'n');
%pthr: for temporal low-pass filter (npixel >= pthr)
pthr = getpar('pthr:',parmat,'n');
%maxsigthr: maximum sigma0
maxsigthr = getpar('maxsigthr:',parmat,'n');

%%% APS filter parameters
%tfiltmethod: aps filter method in time domain (1: time series, 2: p2p stacking)
%xwinsz: for spatial high-pass filter
%ywinsz: for spatial high-pass filter
tfiltmethod = getpar('tfiltmethod:',parmat,'n');
xwinsz = getpar('xwinsz:',parmat,'n');
ywinsz = getpar('ywinsz:',parmat,'n');

%%% reference point coordinate %%%
refx = getpar('refx:',parmat,'n');
refy = getpar('refy:',parmat,'n');

%%% multi-looks
lksx = getpar('lksx:',parmat,'n');
lksy = getpar('lksy:',parmat,'n');

%%% polynomial order for the orbit error fitting
%order = 1: planar; order = 2, quadratic
order = getpar('order:',parmat,'n');

%%% profile parameters
pfaultx0 = getpar('pfaultx0:',parmat,'n');
pfaulty0 = getpar('pfaulty0:',parmat,'n');
pfaultx1 = getpar('pfaultx1:',parmat,'n');
pfaulty1 = getpar('pfaulty1:',parmat,'n');
pfault=[pfaultx0,pfaulty0;pfaultx1,pfaulty1];
swath = getpar('swath:',parmat,'n');
nwin = getpar('nwin:',parmat,'n');

%%%% input interferograms directory %%%%
root = getpar('root:',parmat,'s');
demfile = getpar('demfile:',parmat,'s');
ifg_list = getpar('ifg_list:',parmat,'s');
modelfilename_horz = getpar('modelfilename_horz:',parmat,'s');
obsdir = strcat(root,'obs/');
simdir = strcat(root,'sets/',num2str(simsliprate_horz),'mm/');
modeldir = strcat(root,'fmodels/');
fmodelfile_horz = strcat (modeldir,modelfilename_horz);
demfile = strcat (obsdir,demfile);
ifgfilelist = strcat(obsdir,ifg_list);

%%% output file
outdir = getpar('outdir:',parmat,'s');
ratefile = getpar('ratefile:',parmat,'s');
errorfile = getpar('errorfile:',parmat,'s');
rateproffile = getpar('rateproffile:',parmat,'s');
demproffile = getpar('demproffile:',parmat,'s');
modelproffile = getpar('modelproffile:',parmat,'s');
outdir = strcat (root,outdir);
logfile = strcat(outdir,'log');
ratefile = strcat(outdir,ratefile);
errorfile = strcat(outdir,errorfile);
rateproffile = strcat(outdir,rateproffile);
demproffile = strcat(outdir,demproffile);
modelproffile = strcat(outdir,modelproffile);

%open logfile
logfid = fopen(logfile,'a');

%%% get interferogram list
[ifg_nml,BaseT,ifglist,epochlist] = getnml(ifgfilelist);

%%% make the initial interseismic model
[fmodel_horz]=make_model(fmodelfile_horz,rows,cols);

%%% read dem data %%%%
[dem]=readdem(demfile,rows,cols);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% for vertical slip
if strcmp(rake,'variable')==1
  inisliprate_vert = 1;
  simsliprate_vert = 1;
  modelfilename_vert = getpar('modelfilename_vert:',parmat,'s');
  fmodelfile_vert = strcat (modeldir,modelfilename_vert);
  [fmodel_vert]=make_model(fmodelfile_vert,rows,cols);
else
  inisliprate_vert = 0;
  simsliprate_vert = 0;
  fmodel_vert=zeros(rows,cols);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iterative inverting the slip rate                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for set=1:nsets

  %%% read interferogram
  if realdata==0
    %%% read the synthetic data here
    display(['reading synthetic dataset: ' num2str(set,'%03d')]);
    fprintf(logfid,'reading synthetic dataset: %03d \n',set);
    setdir = strcat(simdir,'set_',num2str(set,'%03d'),'/');
    [ifg]=readifg(setdir,ifg_nml,rows,cols);
  else
    %%% read the real interferograms
    display('reading real data ...');
    fprintf(logfid,'reading real data ...\n');
    [ifg]=readifg(obsdir,ifg_nml,rows,cols);
  end

  %%% find the reference pixel coordinates, run this step only once for each set
  if set==1
    [refx,refy]=refpixel(ifg,refx,refy);
  end

  %%% initialize prevsliprate
  prevsliprate_horz=inisliprate_horz;
  prevsliprate_vert=inisliprate_vert;
  display(['initial slip rate is (h,v): ' num2str(prevsliprate_horz) ' and ' num2str(prevsliprate_vert)]);
  fprintf(logfid,'initial slip rate is (h,v): %-6.2f %-6.2f\n',prevsliprate_horz,prevsliprate_vert);

  %%% iterative estimation for the slip rate
  delta_horz=999;    %initial sliprate changes (horizontal)
  delta_vert=999;    %initial sliprate changes (vertical)
  tol=0.2;           %tolerance for the convergence
  iter=1;            %initial iterative number
  maxiter=10;        %maximam iterative number

  while ((abs(delta_horz)>tol || abs(delta_vert)>tol) && iter<maxiter)
    display(['loop: ' num2str(iter)]);
    fprintf(logfid,'loop: %03d \n',iter);
    
    %% everything is done by loop_fixed
    [fitsliprate_horz,sigma_horz,fitsliprate_vert,sigma_vert,ratemap,errormap,alpha,sigma0]=loop_fixed(ifg,prevsliprate_horz,fmodel_horz,dem,lksx,lksy,ifglist,epochlist,BaseT,order,refx,refy,psize,stackmethod,nsig,pthr,maxsigthr,tfiltmethod,xwinsz,ywinsz,rake,prevsliprate_vert,fmodel_vert);
    
    display(['slip rate is (h,v): ' num2str(fitsliprate_horz) ' +- ' num2str(sigma_horz) ' and ' num2str(fitsliprate_vert) ' +- ' num2str(sigma_vert)]);
    display(['misfit is (sigma0): ' num2str(sigma0) ]);
    fprintf(logfid,'slip rate is (h,v): %6.2f +- %-6.2f %6.2f +- %-6.2f \n',fitsliprate_horz,sigma_horz,fitsliprate_vert,sigma_vert);
    fprintf(logfid,'misfit is (sigma0): %6.2f \n',sigma0);

    delta_horz=fitsliprate_horz-prevsliprate_horz;
    delta_vert=fitsliprate_vert-prevsliprate_vert;
    prevsliprate_horz = fitsliprate_horz;
    prevsliprate_vert = fitsliprate_vert;
    iter=iter+1;
  end
  
  if iter>=maxiter
     fitsliprate_horz=NaN;
     fitsliprate_vert=NaN;
  end

  if realdata==0
    sliprateset_horz(set)=fitsliprate_horz;
    sliprateset_vert(set)=fitsliprate_vert;
    display('---------------------------------------------');
    fprintf(logfid,'---------------------------------------------- \n');
    fprintf(logfid,'final slip rate for synthetic dataset %03d is (h,v): %6.2f +- %-6.2f %6.2f +- %-6.2f\n',set,fitsliprate_horz,sigma_horz,fitsliprate_vert,sigma_vert);
  else
    finalsliprate_horz=fitsliprate_horz;
    finalsliprate_vert=fitsliprate_vert;
    display(['final slip rate is (h,v): ' num2str(fitsliprate_horz) ' +- ' num2str(sigma_horz) ' and ' num2str(fitsliprate_vert) ' +- ' num2str(sigma_vert)]);
    display(['final misfit is (sigma0): ' num2str(sigma0) ]);
    fprintf(logfid,'final slip rate is (h,v): %6.2f +- %-6.2f %6.2f +- %-6.2f\n',fitsliprate_horz,sigma_horz,fitsliprate_vert,sigma_vert);
    fprintf(logfid,'final misfit is (sigma0): %6.2f \n',sigma0);
    % write the output files
    igramfid = fopen(ratefile,'w','l');
    fwrite(igramfid,ratemap','real*4');
    fclose(igramfid);
    igramfid = fopen(errorfile,'w','l');
    fwrite(igramfid,errormap','real*4');
    fclose(igramfid);
    display('calculating profiles...');
    % slip rate adding back the previous model
    ma=mean(alpha);
    [prof_stack,rateprof]=profile(ratemap,pfault,swath,nwin,psize,errormap,ma);
    %plotstack(ratemap,prof_stack,rateprof,1);
    % dem profile
    [prof_stack,demprof]=profile(dem,pfault,swath,nwin,psize);
    % model profile
    fmodel_horz(isnan(ratemap))=NaN;
    [prof_stack,modelprof_horz]=profile(fmodel_horz,pfault,swath,nwin,psize);
    modelprof_horz(:,2)=modelprof_horz(:,2).*finalsliprate_horz;

    % output profile
    proffid = fopen(rateproffile,'w');
    fprintf(proffid,'distance \t rate \t sig_rate \n');
    fprintf(proffid,'%-6.2f %-6.2f %-6.2f \n',rateprof');
    fclose(proffid);
    proffid = fopen(demproffile,'w');
    fprintf(proffid,'distance \t topo \t sig_topo \n');
    fprintf(proffid,'%-6.2f  %-6.2f %-6.2f \n',demprof');
    fclose(proffid);
    proffid = fopen(modelproffile,'w');
    fprintf(proffid,'distance \t model \t sig_model\n');
    fprintf(proffid,'%-6.2f %-6.2f %-6.2f\n',modelprof_horz');
    fclose(proffid);
    break;
  end
end

if realdata==0
  sliprateset_horz(isnan(sliprateset_horz))=[]; %delete non-convergent solutions
  sliprateset_vert(isnan(sliprateset_vert))=[]; %delete non-convergent solutions
  finalsliprate_horz = mean(sliprateset_horz);
  errorsliprate_horz = std(sliprateset_horz);
  finalsliprate_vert = mean(sliprateset_vert);
  errorsliprate_vert = std(sliprateset_vert);
  display(num2str(sliprateset));
  display(['final slip rate is (h,v): ' num2str(finalsliprate_horz) ' +- ' num2str(errorsliprate_horz) ' and ' num2str(finalsliprate_vert) ' +- ' num2str(errorsliprate_vert)]);
  display(['final convergent number is: ' num2str(length(sliprateset_horz))]);
  fprintf(logfid,'final slip rate is (h,v): %6.2f +- %-6.2f  %6.2f +- %-6.2f \n',finalsliprate_horz,errorsliprate_horz,finalsliprate_vert,errorsliprate_vert);
  fprintf(logfid,'final convergent number is: %4d \n',length(sliprateset_horz));
end

fclose(logfid);
display('finished...');
