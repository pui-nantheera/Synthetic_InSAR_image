%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  atmfilter.m: atm filter in time and space domain using PS method  %
%%  author: Hua Wang                                                  %
%%  date: 02/02/2008                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ifg_aps,ifg_stack,std_stack]=atmfilter(ifg,ifglist,epochlist,BaseT,vcm_t,pthr,nsig,tfiltmethod,xwinsz,ywinsz)

[nifgs]=size(ifg,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temporal low-pass filter                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ifg_apst: atmospheric phase screen (APS) in time domain

if tfiltmethod==1
  % temporal low-pass filter using time series analysis method %
  r2=5;
  [ifg_ts,ifg_fit,ifg_apst,ifg_var,ts_rough]=tsinvlap(ifg,ifglist,epochlist,vcm_t,pthr,r2);
  %[ifg_ts]=tsinv(ifg,ifglist,epochlist,vcm_t,pthr);
else
  % temporal low-pass filter using stacking method %
  [ifg_stack,std_stack]=stack_p2p(ifg,BaseT,vcm_t,pthr,nsig);
  for i=1:nifgs
    ifg_apst(:,:,i) = ifg(:,:,i)-ifg_stack*BaseT(i);
  end
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spatial high-pass filter                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ifg_aps: atmospheric phase screen (APS)
for i=1:nifgs
  ifg_aps(:,:,i) = meanfilter(ifg_apst(:,:,i),xwinsz,ywinsz);
end
