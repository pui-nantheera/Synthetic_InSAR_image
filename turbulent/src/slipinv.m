function [s,stdx,tiltmap,resmap,sigma0]=slipinv(ratemap,errormap,fmodel_horz,vcm_s,refx,refy,lksx,lksy,psize,rake,fmodel_vert)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%re-calculate the full variance-covariance matrix of the rate map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ratemap_lowres] = looks(ratemap,lksx,lksy);
[errormap_lowres] = looks(errormap,lksx,lksy);
[fmodel_h_lowres] = looks(fmodel_horz,lksx,lksy);
[rows_lowres,cols_lowres]=size(ratemap_lowres);
n_lowres=rows_lowres*cols_lowres;
refx_lowres=floor(refx/lksx);
refy_lowres=floor(refy/lksy);

sig = reshape(errormap_lowres',n_lowres,1);
vcm_tmp = sig*sig';
vcm_s = vcm_tmp.*vcm_s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%invert the model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datav=reshape(ratemap_lowres',n_lowres,1);
modelv=reshape(fmodel_h_lowres',n_lowres,1);
if strcmp(rake,'variable')==1
  [fmodel_v_lowres] = looks(fmodel_vert,lksx,lksy);
  modelv_v=reshape(fmodel_v_lowres',n_lowres,1);
  modelv=[modelv modelv_v];
end

[yy,xx]=meshgrid(1:rows_lowres,1:cols_lowres);
xxv=reshape(xx,n_lowres,1);
yyv=reshape(yy,n_lowres,1);
xxv=(xxv-refx_lowres)*psize*lksx;
yyv=(yyv-refy_lowres)*psize*lksy;
tiltsetup=[xxv yyv ones(n_lowres,1)];
B=[modelv tiltsetup];

B(isnan(datav),:)=[];
vcm_s(isnan(datav),:)=[];
vcm_s(:,isnan(datav))=[];
datav(isnan(datav))=[];

%%% pcg method
%the same results obtained by lscov method
%tol=1e-6;iter=1000;
%y=pcg(vcm_s,datav,tol,iter);
%for n=1:4
%  zn(:,:,n)=[pcg(vcm_s,B(:,n),tol,iter)];
%end
%z=[zn(:,:,1) zn(:,:,2) zn(:,:,3) zn(:,:,4)];
%s=pcg(B'*z,B'*y,tol,iter);
%s=s';

%%%calculate the slip rate using least squares
[s,stdx,mse]=lscov(B,datav,vcm_s);
sigma0=sqrt(mse);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%calculate the tilts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npar=length(s);
tiltpar=s(npar-2:npar);
[rows,cols]=size(ratemap);
[yy,xx]=meshgrid(1:rows,1:cols);
xxv=reshape(xx,rows*cols,1);
yyv=reshape(yy,rows*cols,1);
xxv=(xxv-refx)*psize;
yyv=(yyv-refy)*psize;
tiltsetup=[xxv yyv ones(rows*cols,1)];
tiltmap = tiltsetup*tiltpar;
tiltmap = (reshape(tiltmap,cols,rows))';
tiltmap(isnan(ratemap))=NaN;

resmap = s(1)*fmodel_horz+tiltmap-ratemap;
if strcmp(rake,'variable')==1
  resmap = resmap+s(2)*fmodel_vert;
end
