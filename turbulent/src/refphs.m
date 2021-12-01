%%% by Hua Wang on 4 Mar 2008

function[refphs]=refphs(ifg,refx,refy)

lksx=3;
lksy=3;
thr = floor(lksx*lksy/2);

ifg(isnan(ifg))=0;
patch=ifg(refy-1:refy+1,refx-1:refx+1);
mask=ones(lksy,lksx);
mask(patch==0)=0;
n=sum(sum(mask));
refphs=sum(sum(patch))/n;

if n<thr
  disp('Note: the reference pixel is not in high coherent area!'); 
end
