function [ifg_orb,tiltpar,stdtiltpar] = orbfit(ifg,lksx,lksy,refx,refy,psize,order)

[rows,cols,nifgs]=size(ifg);
%% multilook interferogram
for i=1:nifgs
  [ifg_lowres(:,:,i)]=looks(ifg(:,:,i),lksx,lksy);
end 
[rows_lowres,cols_lowres,nifgs]=size(ifg_lowres);

%x,y coordinates
[yy,xx]=meshgrid(1:rows_lowres,1:cols_lowres);
for i=1:nifgs
  %observations
  ifgv = reshape(ifg_lowres(:,:,i)',rows_lowres*cols_lowres,1);
  obsv = ifgv;
  obsv(isnan(ifgv))=[];

  %x,y coordinates
  xxv=reshape(xx,rows_lowres*cols_lowres,1);
  yyv=reshape(yy,rows_lowres*cols_lowres,1);
  xxv(isnan(ifgv))=[];
  yyv(isnan(ifgv))=[];
  %subtract the reference coordinates
  xxv = (xxv-floor(refx/lksx))*psize*lksx;
  yyv = (yyv-floor(refy/lksy))*psize*lksy;
  %coefficient matrix
  if order==1
    B=[xxv yyv ones(length(xxv),1)];
  else
  %order == 2
    xxv2 = xxv.*xxv;
    yyv2 = yyv.*yyv;
    xyv  = xxv.*yyv;
    B=[xxv2 yyv2 xyv xxv yyv ones(length(xxv),1)];
  end

  %least-squares inversion
  %using the unit weight
  [tiltpar(:,i),stdtiltpar(:,i)]=lscov(B,obsv,diag(ones(length(obsv),1)));

end

%%%output the original full resolution interferogram
[yy,xx]=meshgrid(1:rows,1:cols);
%x,y coordinates
xxv=reshape(xx,rows*cols,1);
yyv=reshape(yy,rows*cols,1);
%subtract the reference coordinates
xxv = (xxv-refx)*psize;
yyv = (yyv-refy)*psize;
%coefficient matrix
if order==1
  B=[xxv yyv ones(length(xxv),1)];
else
%order == 2
  xxv2 = xxv.*xxv;
  yyv2 = yyv.*yyv;
  xyv  = xxv.*yyv;
  B=[xxv2 yyv2 xyv xxv yyv ones(length(xxv),1)];
end

for i=1:nifgs
  fit = B*tiltpar(:,i);
  ifg_orb(:,:,i)=(reshape(fit,cols,rows))';
end
ifg_orb(isnan(ifg))=NaN;
