%%pcmc_epoch.m
%synthesize using multi-interferogram method (epoch by epoch)

function pcmc_epoch(simdir,ifg_model,ifg_nml,ifglist,maxvar,alpha,covmodel_type,nsets,psize,lksx,lksy,refx,refy,tiltpar,order)

[rows,cols,nifgs]=size(ifg_model);
nepoch = max(max(ifglist));

rows_lowres=floor(rows/lksy);
cols_lowres=floor(cols/lksx);
psizex=psize*lksx;
psizey=psize*lksy;

%calculate the mean values of orbit/atmospheric parameters
maxvar_mean = mean(maxvar)/2;
alpha_mean = mean(alpha);
tiltpar_mean = mean(tiltpar');
stdtiltpar = std(tiltpar');

%make grid
[yy,xx]=meshgrid(1:rows,1:cols);
[yy_lowres,xx_lowres]=meshgrid(1:lksy:lksy*rows_lowres,1:lksx:lksx*cols_lowres);

%synthesizing
for i=1:nepoch

  display(strcat('synthesizing epoch: ',num2str(i),'/',num2str(nepoch)));

  % using epoch dependent parameters
  %[orb_pets] = pcmc_orb(rows,cols,nsets,psize,refx,refy,tiltpar(:,i)',stdtiltpar,order);
  %[atm_pets] = pcmc_atm(rows_lowres,cols_lowres,maxvar(i),alpha(i),covmodel_type,nsets,psizex,psizey);

  % using the same parameters for all epochs
  [orb_pets] = pcmc_orb(rows,cols,nsets,psize,refx,refy,tiltpar_mean,stdtiltpar,order);
  [atm_pets] = pcmc_atm(rows_lowres,cols_lowres,maxvar_mean,alpha_mean,covmodel_type,nsets,psizex,psizey);

  for iset = 1:nsets
    %% interpolation for each interferogram in each set for atm_pets
    zz_lowres=atm_pets(:,:,iset);
    atm_intp = interp2(xx_lowres',yy_lowres',zz_lowres,xx',yy','bilinear');
   
    %% add peturb phase together
    pha_pets = orb_pets(:,:,iset) + atm_intp;

    %% write out data perturb phase to save memory (make program slower, any other idea???)
    setdir = strcat(simdir,'set_',num2str(iset,'%03d'),'/');
    outname = char(strcat(setdir,'epoch_',num2str(i,'%03d')));
    fid = fopen(outname,'w');
    fwrite(fid,pha_pets','real*4');
    fclose(fid);
  end

  clear orb_pets;
  clear atm_pets;
end

for iset=1:nsets

  if mod(iset,10)==0
    display(strcat('synthesizing set: ',num2str(iset),'/',num2str(nsets)));
  end

  setdir = strcat(simdir,'set_',num2str(iset,'%03d'),'/');

  %% read data epoch by epoch
  for i=1:nepoch
    iname = char(strcat(setdir,'epoch_',num2str(i,'%03d')));
    fin = fopen(iname,'r');
    temp = fread(fin, [cols,rows], 'real*4');
    fclose(fin);
    pha_pets(:,:,i) = temp';
  end

  %generate interferograms
  for i=1:nifgs
    %master/slave image epoch, using it to determine the cofficient position
    im = ifglist(i,1);
    is = ifglist(i,2);

    %add back the model/orbit/atmosphere together
    ifg_pets = ifg_model(:,:,i) - pha_pets(:,:,im) + pha_pets(:,:,is);

    %write output file
    outname = char(strcat(setdir,ifg_nml(i)));
    fid = fopen(outname,'w');
    fwrite(fid,ifg_pets','real*4');
    fclose(fid);
  end
end
