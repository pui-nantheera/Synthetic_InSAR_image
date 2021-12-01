%%% by Hua Wang on 4 Mar 2008

function[ifg_lowres]=looks(ifg,lksx,lksy)

if (lksx==1 && lksy==1)
  ifg_lowres=ifg;
else
  [rows,cols]=size(ifg);
  rows_lowres=floor(rows/lksy);
  cols_lowres=floor(cols/lksx);
  ifg_lowres=zeros(rows_lowres,cols_lowres);

  thr = floor(lksx*lksy/2);

  ifg(isnan(ifg))=0;
  for r=1:rows_lowres 
    for c=1:cols_lowres
      patch=ifg((r-1)*lksy+1:r*lksy,(c-1)*lksx+1:c*lksx);
      mask=ones(lksy,lksx);
      mask(patch==0)=0;
      n=sum(sum(mask));
      if n<thr
         ifg_lowres(r,c)=NaN;
      else
         ifg_lowres(r,c)=sum(sum(patch))/n;
      end
    end
  end
end
