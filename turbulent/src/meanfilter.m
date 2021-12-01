%%% by Hua Wang on 4 Mar 2008

function[ifg_filt]=meanfilter(ifg,xwsz,ywsz)

%xwsz,ywsz: x and y window size

[rows,cols]=size(ifg);

thr = floor(xwsz*ywsz/2);

xwszd2 = floor(xwsz/2);
ywszd2 = floor(ywsz/2);

ifg_filt=ifg;  %paded in the boundary

ifg(isnan(ifg))=0;

for r=ywszd2+1:rows-ywszd2
  for c=xwszd2+1:cols-xwszd2
    patch=ifg(r-ywszd2:r+ywszd2,c-xwszd2:c+xwszd2);
    mask=ones(ywsz,xwsz);
    mask(patch==0)=0;
    n=sum(sum(mask));
    if(n<thr)
      ifg_filt(r,c)=NaN;
    else
      ifg_filt(r,c)=sum(sum(patch))/n;
    end
  end
end
