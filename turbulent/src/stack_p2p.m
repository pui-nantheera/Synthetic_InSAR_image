%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  stack_p2p.m: stacking pixel by pixel      %
%%  author: Hua Wang                          %
%%  date: 02/02/2008                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ifg_stack,std_stack]=stack_p2p(ifg,BaseT,vcm,pthr,nsig,maxsigthr)

[rows,cols,nifgs]=size(ifg);

ifg = permute(ifg,[3,1,2]);
std_stack=NaN(rows,cols);
ifg_stack=NaN(rows,cols);

for i=1:rows
  for j=1:cols
    %initialize
    ifgv = ifg(:,i,j);
    mask=ones(nifgs,1);
    mask(isnan(ifgv))=0;
    m=sum(mask);
    if m>=pthr
      vcm_tmp = vcm;
      B = BaseT;
      %delete NaN pixels
      B(isnan(ifgv))=[];
      vcm_tmp(isnan(ifgv),:)=[];
      vcm_tmp(:,isnan(ifgv))=[];
      ifgv(isnan(ifgv))=[];

      %calculate slip rate by iterative least-squares
      while (m>=pthr)
        %%least squares stacking
        %[x,stdx,mse]=lscov(B,ifgv,vcm_tmp);

        P=pinv(vcm_tmp);          %calculates inverse vcm.
        Nbb=B'*P*B;             
        Nbb=1/Nbb;
        W=B'*P*ifgv;            
        x=Nbb*W;                  %calculates rate on pixel
      
        %calculate residuals
        v = B*x-ifgv;
        mse = v'*P*v/(length(ifgv)-1);
        sig0 = sqrt(mse);
        stdx = sig0*sqrt(Nbb);

        %delete the maximum outliers 
        [maxv,maxi] = max(abs(v));
        if (maxv>nsig*sig0)
          m=m-1;
          if m>=pthr
            B(maxi)=[];
            vcm_tmp(maxi,:)=[];
            vcm_tmp(:,maxi)=[];
            ifgv(maxi)=[];
          end
        else
          ifg_stack(i,j)=x;
          std_stack(i,j)=stdx;
          break;
        end
      end
    end
  end
end

%mask the data whose sigma is larger than the maxmum sigma threshold
ifg_stack(std_stack>maxsigthr) = NaN;
std_stack(std_stack>maxsigthr) = NaN;
