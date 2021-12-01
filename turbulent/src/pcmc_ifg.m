%%pcmc_ifg.m
%%synthesize using interferogram dependent method

function pcmc_ifg(simdir,ifg_model,ifg_nml,maxvar,alpha,covmodel_type,nsets,psize,lksx,lksy,refx,refy,tiltpar,stdtiltpar,order)

[rows,cols,nifgs]=size(ifg_model);

rows_lowres=floor(rows/lksy);
cols_lowres=floor(cols/lksx);
psizex=psize*lksx;
psizey=psize*lksy;

%make grid
[yy,xx]=meshgrid(1:rows,1:cols);
[yy_lowres,xx_lowres]=meshgrid(1:lksy:lksy*rows_lowres,1:lksx:lksx*cols_lowres);   

for i=1:nifgs

  display(strcat('synthesizing interferogram: ',num2str(i),'/',num2str(nifgs)));

  %% 1. for orbit error
  [orb_pets] = pcmc_orb(rows,cols,nsets,psize,refx,refy,tiltpar(:,i),stdtiltpar(:,i),order);
  
  %% 2. for atmospheric delay
  [atm_pets] = pcmc_atm(rows_lowres,cols_lowres,maxvar(i),alpha(i),covmodel_type,nsets,psizex,psizey);

  for iset=1:nsets
    %% 3. interpolation for each interferogram in each set for atm_pets
    zz_lowres=atm_pets(:,:,iset);
    atm_intp = interp2(xx_lowres',yy_lowres',zz_lowres,xx',yy','bilinear');

    % 4. add back the model/orbit/atmosphere together
    pha_pets = ifg_model(:,:,i) + orb_pets(:,:,iset) + atm_intp;

    %write output file
    setdir = strcat(simdir,'set_',num2str(iset,'%03d'),'/');
    outname = char(strcat(setdir,ifg_nml(i)));
    fin = fopen(outname,'w');
    fwrite(fin,pha_pets','real*4');
    fclose(fin);
  end
end
