%%pcmc_orb.m

function [orb_pets] = pcmc_orb(rows,cols,N,psize,refx,refy,tiltpar,stdtiltpar,order)

n = rows*cols;
ncoef=(order+1)*(order+2)/2-1;

%1. create matrix of N gaussian noise vectors length n (on it's side!)
for i=1:ncoef
  Z = randn(1,N);
  simtiltpar(i,:)=tiltpar(i)+Z.*stdtiltpar(i);
end

%2. design matrix
%make grid
[yy,xx]=meshgrid(1:rows,1:cols);
xxv=reshape(xx,n,1);
yyv=reshape(yy,n,1);
%subtract the reference coordinates
xxv = (xxv-refx)*psize;
yyv = (yyv-refy)*psize;

if order==1
  %B=[xxv yyv ones(length(xxv),1)];
  B=[xxv yyv];
else
  %order == 2
  xxv2 = xxv.*xxv;
  yyv2 = yyv.*yyv;
  xyv  = xxv.*yyv;
  %B=[xxv2 yyv2 xyv xxv yyv ones(length(xxv),1)];
  B=[xxv2 yyv2 xyv xxv yyv];
end

%3. peturbed phase
for i = 1:N
  X = B*simtiltpar(:,i);
  pha_pets = reshape(X,cols,rows);
  orb_pets(:,:,i) = pha_pets';
end
