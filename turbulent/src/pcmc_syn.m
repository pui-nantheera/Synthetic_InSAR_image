%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  pcmc_syn.m: main program for synthesization               %
%%  date: 13/04/2008                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=pcmc_syn(parfile)

%read parameter file
[parmat] = readparfile(parfile);

%%%% interferogram size %%%%%%%%
rows  = getpar('rows:',parmat,'n');
cols  = getpar('cols:',parmat,'n');
psize = getpar('psize:',parmat,'n');

%%%% inisliprate slip rate for the observed interferogram in mm/yr
inisliprate_horz = getpar('inisliprate_horz:',parmat,'n');
%%%% synthetic slip rate in mm/yr
simsliprate_horz = getpar('simsliprate_horz:',parmat,'n');

%%% number of sets of data (multiple sets required for Monte Carlo error estimates %%%
nsets = getpar('nsets:',parmat,'n');

%%% stacking method
stackmethod = getpar('stackmethod:',parmat,'n');

%%% reference point coordinate %%%
refx = getpar('refx:',parmat,'n');
refy = getpar('refy:',parmat,'n');

%%% multi-looks
lksx = getpar('lksx:',parmat,'n');
lksy = getpar('lksy:',parmat,'n');

%%% polynomial order for the orbit error fitting
%order = 1: planar; order = 2, quadratic
order = getpar('order:',parmat,'n');


%%%% synthetic method
%simmethod = 2; %1: interferogram dependent method; 2: multi-interferogram method (recommended)
simmethod = getpar('simmethod:',parmat,'n');

%%%% input interferograms directory %%%%
root = getpar('root:',parmat,'s');
fg_list = getpar('ifg_list:',parmat,'s');
obsdir = strcat(root,'obs/');
simdir = strcat(root,'sets/');
simdir = strcat(simdir,num2str(simsliprate),'mm/');
modeldir = strcat(root,'fmodels/');
modelfilename_horz = getpar('modelfilename_horz:',parmat,'s');
fmodelfile_horz = strcat (modeldir,modelfilename_horz);
ifgfilelist = strcat(obsdir,ifg_list);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate the vcm of atmospheric delay and orbit errors from the real data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% get interferogram list
[ifg_nml,BaseT,ifglist,epochlist] = getnml(ifgfilelist);

%%% make the initial interseismic model
[fmodel]=make_model(fmodelfile_horz,rows,cols);

%%% read the real interferograms
[ifg]=readifg(obsdir,ifg_nml,rows,cols);

%%% find the reference pixel coordinates
[refx,refy]=refpixel(ifg,refx,refy);

%%% make the initial interseismic model
nifgs=length(BaseT);
[fmodel]=make_model(fmodelfile,rows,cols);
for i=1:nifgs
  ifg_model(:,:,i)=inisliprate*fmodel*BaseT(i);
end
ifg_model(isnan(ifg))=NaN;
ifg_nomodel = ifg-ifg_model;

%%% polynomial fitting the orbit error
[ifg_orb,tiltpar] = orbcorrect(ifg_nomodel,ifglist,order,lksx,lksy,refx,refy,psize,stackmethod);
ifg_flat=ifg_nomodel-ifg_orb;

%%% calculate covariance matrix for each interferogram
[vcm_t,vcm_s,maxvar,alpha] = make_vcm(ifg_flat,ifglist,psize,lksx,lksy);

%%% clear some local variables
clear ifg_nomodel;
clear ifg_orb;
clear ifg_flat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% synthesize the interferograms here                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('synthesizing interferogram...');
ifg_model = ifg_model .* (simsliprate/inisliprate);
if simmethod==1
  pcmc_ifg(simdir,ifg_model,ifg_nml,maxvar,alpha,0,nsets,psize,lksx,lksy,refx,refy,tiltpar,stdtiltpar,order);
else
  pcmc_epoch(simdir,ifg_model,ifg_nml,ifglist,maxvar,alpha,0,nsets,psize,lksx,lksy,refx,refy,tiltpar,order);
end
