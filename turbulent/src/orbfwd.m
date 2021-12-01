function [ifg_orb]=orbfwd(ifglist,order,rows,cols,refx,refy,psizex,psizey,stackmethod,orbparams)

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

%x,y coordinates
[yy,xx]=meshgrid(1:rows,1:cols);
xxv=reshape(xx,rows*cols,1);
yyv=reshape(yy,rows*cols,1);
%subtract the reference coordinates
xxv = (xxv-refx)*psizex;
yyv = (yyv-refy)*psizey;
clear xx;
clear yy;

if order==1
  B=[xxv yyv ones(length(xxv),1)];
else
%order == 2
  xxv2 = xxv.*xxv;
  yyv2 = yyv.*yyv;
  xyv  = xxv.*yyv;
  B=[xxv2  yyv2  xyv  xxv  yyv ones(length(xxv),1)];
end
clear xxv;
clear yyv;

for i=1:nifgs
  %multi-interferogram method matrix
  if stackmethod==3
    %master/slave image epoch, using it to determine the cofficient position
    im = ifglist(i,1);
    is = ifglist(i,2);
    %coefficients for the master/slave image
    jbm=(im-1)*ncoef+1;
    jbs=(is-1)*ncoef+1;
    parmm = orbparams(jbm:jbm+ncoef-1);
    parms = orbparams(jbs:jbs+ncoef-1);
    parm = parms - parmm;
  %ifg by ifg method matrix
  else
    jb1 = (i-1)*ncoef+1;
    jb2 = i*ncoef;
    parm = orbparams(jb1:jb2);
  end

  %coefficients for the offsets
  joff = npolypar+i;
  parm = [parm;orbparams(joff)];

  fullorb = B*parm;
  ifg_orb(:,:,i)=reshape(fullorb,cols,rows)';
end
