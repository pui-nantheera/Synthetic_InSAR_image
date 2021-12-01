function [ifg_dem] = atmfit(ifg,dem,lksx,lksy)

[rows,cols,nifgs]=size(ifg);

%% multilook interferogram
for i=1:nifgs
  [ifg_lowres(:,:,i)]=looks(ifg(:,:,i),lksx,lksy);
end 
[dem_lowres]=looks(dem,lksx,lksy);
[rows_lowres,cols_lowres,nifgs]=size(ifg_lowres);

for i=1:nifgs
  %observations
  ifgv = reshape(ifg_lowres(:,:,i)',rows_lowres*cols_lowres,1);
  obsv = ifgv;
  obsv(isnan(ifgv))=[];

  %demdata
  demv = reshape(dem_lowres',rows_lowres*cols_lowres,1);
  demv(isnan(ifgv))=[];

  %coefficient matrix
  B=[demv ones(length(demv),1)];

  %least-squares inversion
  %using the unit weight
  invnbb = inv(B'*B);
  tiltpar(:,i) = invnbb*B'*obsv;
end

%%%output the original full resolution interferogram
demv = reshape(dem',rows*cols,1);

for i=1:nifgs
  B=[demv ones(length(demv),1)];
  fit = B*tiltpar(:,i);
  ifg_dem(:,:,i)=(reshape(fit,cols,rows))';
end
ifg_dem(isnan(ifg))=NaN;
