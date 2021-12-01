function [sliprate_horz,sigma_horz,sliprate_vert,sigma_vert,ratemap,errormap,alpha,sigma0]=loop_fixed(ifg,prevsliprate_horz,fmodel_horz,dem,lksx,lksy,ifglist,epochlist,BaseT,order,refx,refy,psize,stackmethod,nsig,pthr,maxsigthr,tfiltmethod,xwinsz,ywinsz,rake,prevsliprate_vert,fmodel_vert);

%%%dimension of ifg
nifgs=length(ifglist);

% remove reference value from all igrams (removed this step on 30/04/2008)
% !! useless, the phase unwrapping offset can be considered during orbit correction

%load mask
%mask(mask==0)=NaN;
%for i=1:nifgs
%  ifg(:,:,i)=ifg(:,:,i).*mask;
%end

%remove initial model
disp('removing the initial model...');
for i=1:nifgs
  ifg_model(:,:,i)=(prevsliprate_horz*fmodel_horz+prevsliprate_vert*fmodel_vert)*BaseT(i);
end
ifg=ifg-ifg_model;

%correct orbit errors
disp('correcting orbit errors...');
[ifg_orb,tiltpar] = orbcorrect(ifg,ifglist,order,lksx,lksy,refx,refy,psize,stackmethod);
ifg=ifg-ifg_orb;
clear ifg_orb;

%%% calculate covariance matrix for each interferogram
[vcm_t,vcm_s,maxvar,alpha] = make_vcm(ifg,ifglist,psize,lksx,lksy);
%[vcm_t,vcm_s] = make_vcm_correct(ifg,ifglist,epochlist,psize,lksx,lksy);

%correct atmospheric delay errors
%disp('correcting atmospheric delay errors...');
%[ifg_aps,ifg_stack,std_stack]=atmfilter(ifg,ifglist,epochlist,BaseT,vcm_t,pthr,nsig,tfiltmethod,xwinsz,ywinsz);
%ifg=ifg-ifg_aps;

%%% correcting DEM derived atmospheric delay errors
[ifg_atm] = atmcorrect(ifg,ifglist,dem,lksx,lksy,refx,refy,stackmethod);
ifg=ifg-ifg_atm;
clear ifg_atm;

%%% add back the initial model
ifg=ifg+ifg_model;
clear ifg_model;

%%%Using Hua Wang's program to stack pixel by pixel
disp('stacking (it may take a few minutes) ...');
if stackmethod==1
  [ratemap]=stack_cls(ifg,BaseT);
else
  [ratemap,errormap]=stack_p2p(ifg,BaseT,vcm_t,pthr,nsig,maxsigthr);
end

%%%mask some patch
%load mask
%mask(mask==0)=NaN;
%ratemap=ratemap.*mask;
%errormap=errormap.*mask;

%%%estimates slip rate on fault from rate map.
disp('inverting slip rate...');
%[s,stdx,tiltmap,resmap,sigma0]=slipinv(ratemap,errormap,fmodel,vcm_s,refx,refy,lksx,lksy,psize,rake,fmodel_vert);
 [s,stdx,tiltmap,atmmap,resmap,sigma0]=slipinvdem(ratemap,errormap,fmodel_horz,dem,vcm_s,refx,refy,lksx,lksy,psize,rake,fmodel_vert);
sliprate_horz=s(1);
sigma_horz=stdx(1);
if strcmp(rake,'variable')==1
  sliprate_vert=s(2);
  sigma_vert=stdx(2);
else
  sliprate_vert=0;
  sigma_vert=0;
end
%%%ratemap without orbit and atm error
ratemap = sliprate_horz*fmodel_horz + sliprate_vert*fmodel_vert - resmap;
