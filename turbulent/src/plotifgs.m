%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  plotifgs.m: plot multiple interferograms  %
%%  date: 02/02/2008                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=plotifgs(ifg,m,n,ifg_nml)

nifgs=size(ifg,3);

%% plot the interferogram
figure

for i=1:nifgs
  if nargin>3
    strpair=char(ifg_nml(i)); 
    title(strpair(5:17));
  end
  subplot(m,n,i);
  imagesc(ifg(:,:,i));  
  colormap(jet);
  axis equal;
  axis off;
  box on;
  colorbar;

  %Now set the alpha map for the NaN region
  z = ifg(:,:,i);
  z(~isnan(ifg(:,:,i))) = 1;
  z(isnan(ifg(:,:,i))) = 0;
  alpha(z);
  set(gca, 'color', [0 0 0]);
end
%colorbar('location','southoutside');
