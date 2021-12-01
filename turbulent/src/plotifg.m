%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  plotifg.m: plot one interferogram         %
%%  author: Hua Wang                          %
%%  date: 02/02/2008                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=plotifg(ifg,clim,strtitle)

[rows,cols]=size(ifg);

if nargin<2
  clim=2;
end

%% plot the interferogram
%figure(iplot)
figure

%calculate std and mean
ifgv = reshape(ifg',rows*cols,1);
ifgv(isnan(ifgv))=[];
std_ifg=std(ifgv);
mean_ifg=mean(ifgv);

%display(num2str(mean_ifg));
%display(num2str(std_ifg));

if clim==1
  clims = [mean_ifg-2*std_ifg mean_ifg+2*std_ifg];
else
  clims = [min(min(ifg)) max(max(ifg))];
end

if nargin==3
  title(strtitle);
end
imagesc(ifg,clims);
colormap(jet);
colorbar;
axis equal;
axis off;
box on;

%Now set the alpha map for the NaN region
z = ifg(:,:);
z(~isnan(ifg)) = 1;
z(isnan(ifg)) = 0;
alpha(z);
set(gca, 'color', [0 0 0]);

%filename = strcat(num2str(iplot),'.jpg');
%print('-djpeg', filename);
