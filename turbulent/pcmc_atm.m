%%pcmc.m  Peturb using Correlated Monte Carlo
%%
%%Matlab script to perturb interferogram files using a given cov. function
%% Type 2 'ebessel' recommended
%%
%% tjw 2004

%covmodel_type = 0 %0=exponential; 1=expcos; 2=ebessel
function [atm_pets] = pcmc_atm(rows,cols,maxvar,alpha,covmodel_type,N,psizex,psizey)

n = rows*cols;
%%calculate distance between points in a matrix
[yy,xx]=meshgrid(1:rows,1:cols);
xxv=reshape(xx,n,1);
yyv=reshape(yy,n,1);
xxv=xxv.*psizex;
yyv=yyv.*psizey;
dx=repmat(xxv,1,n)-repmat(xxv',n,1);
dy=repmat(yyv,1,n)-repmat(yyv',n,1);
rgrid=sqrt(dx.^2+dy.^2);

clear xx;
clear yy;
clear xxv;
clear yyv;
clear dx;
clear dy;

if (covmodel_type==0)
  %%Use exp function to calculate vcm
  vcm = maxvar*exp(-alpha*rgrid);
elseif (covmodel_type==1)
  %%Use expcos function to calculate vcm
  vcm = maxvar*exp(-alpha*rgrid).*cos(beta*rgrid);
elseif (covmodel_type==2)
  vcm = ebessel(rgrid,maxvar,eb_r,eb_w);
end

%%Calculate correlated noise using Cholesky Decomposition
%1. create matrix of N gaussian noise vectors length n (on it's side!) 
Z = randn(N,n);
%2. chol decomp on vcm
V = chol(vcm); %upper triangular part of cholesky [V'V=vcm]
%3. Create matrix X containing N correlated noisevecotrs length n (on it's side!)
X = Z * V;
X = X'; %transpose to make it N vectors of length n the right way up

%4. write out peturbed files
for i = 1:N
  pha_pets = reshape(X(:,i),cols,rows);
  atm_pets(:,:,i) = pha_pets';
end
