%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  stack_cls.m: interferogram stacking       %
%%  author: Hua Wang                          %
%%  date: 02/02/2008                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[ifg_stack]=stack_cls(ifg,BaseT)

[rows,cols,nifgs]=size(ifg);
%%% stacking the public pixels in all interferograms
ifg_stack=zeros(rows,cols);
sumt = sum(BaseT);
for i=1:nifgs
  ifg_stack=ifg_stack+ifg(:,:,i);
end
ifg_stack = ifg_stack/sumt;
%end
