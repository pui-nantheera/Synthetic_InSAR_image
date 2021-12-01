function [B]=orbdesign(ifglist,order,rows,cols,refx,refy,psizex,psizey,stackmethod)

nifgs=length(ifglist);

%cofficients number for each epoch
%order=1: x1, y1
%order=2: x1^2, y1^2, x1y1, x1, y1
ncoef = (order+1)*(order+2)/2-1;

%epoch number
nepoch = max(max(ifglist));

%total parameters number (epoch parameters + offsets)
if stackmethod==3
  npolypar = nepoch*ncoef; %multi-interferogram method
else
  npolypar = nifgs*ncoef; %ifg by ifg method
end
npar = npolypar + nifgs;

%x,y coordinates
[yy,xx]=meshgrid(1:rows,1:cols);
xxv=reshape(xx,rows*cols,1);
yyv=reshape(yy,rows*cols,1);
%subtract the reference coordinates
xxv = (xxv-refx)*psizex;
yyv = (yyv-refy)*psizey;
clear xx;
clear yy;

B = zeros(length(xxv)*nifgs,npar);
for i=1:nifgs
  %start/end line number for each interferogram in the coefficient matrix
  ib1 = (i-1)*length(xxv)+1;
  ib2 = i*length(xxv);

  %multi-interferogram method matrix
  if stackmethod==3
    %master/slave image epoch, using it to determine the cofficient position
    im = ifglist(i,1);
    is = ifglist(i,2);
    %coefficients for the master/slave image
    jbm=(im-1)*ncoef+1;
    jbs=(is-1)*ncoef+1;
    if order==1
      B(ib1:ib2,jbm:jbm+ncoef-1)=[-xxv -yyv];
      B(ib1:ib2,jbs:jbs+ncoef-1)=[ xxv  yyv];
    else
    %order == 2
      xxv2 = xxv.*xxv;
      yyv2 = yyv.*yyv;
      xyv  = xxv.*yyv;
      B(ib1:ib2,jbm:jbm+ncoef-1)=[-xxv2 -yyv2 -xyv -xxv -yyv];
      B(ib1:ib2,jbs:jbs+ncoef-1)=[ xxv2  yyv2  xyv  xxv  yyv];
    end

  %ifg by ifg method matrix
  else
    jb1 = (i-1)*ncoef+1;
    jb2 = i*ncoef;
    if order==1
      B(ib1:ib2,jb1:jb2)=[xxv yyv];
    else
    %order == 2
      xxv2 = xxv.*xxv;
      yyv2 = yyv.*yyv;
      xyv  = xxv.*yyv;
      B(ib1:ib2,jb1:jb2)=[xxv2 yyv2 xyv xxv yyv];
    end
  end

  %coefficients for the offsets
  joff = npolypar+i;
  B(ib1:ib2,joff)=ones(length(xxv),1);
end
