function [refx,refy] = refpixel(ifg,refx,refy)

[rows,cols,nifgs]=size(ifg);

% find reference pixel
%choose first coherent pixel in list as reference pixel
if (refx==0) && (refy==0)
  for r=1:rows
    for c=1:cols
      ncoh = 0;
      for i=1:nifgs
        if isnan(ifg(r,c,i))==1
          break;
        else
          ncoh = ncoh+1;
        end
      end
      if ncoh==nifgs
        refx=c;
        refy=r;
        break;
      end
    end
    if (refx~=0) && (refy~=0)
      break;
    end
  end
end
