%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  plotstack.m: plot stacked interferogram   %
%%  author: Hua Wang                          %
%%  date: 02/02/2008                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=plotstack(ifg,prof,sumprof,clim)

[rows,cols]=size(ifg);

figure

%calculate std and mean
ifgv = reshape(ifg',rows*cols,1);
obsv = ifgv;
obsv(isnan(ifgv))=[];
std_ifg=std(obsv);
mean_ifg=mean(obsv);

%display(num2str(mean_ifg));
%display(num2str(std_ifg));

if clim==1
  clims = [mean_ifg-2*std_ifg mean_ifg+2*std_ifg];
else
  clims = [min(min(ifg)) max(max(ifg))];
end

subplot(1,2,1);
imagesc(ifg,clims);  
colormap(jet);
colorbar('Location','SouthOutSide');
axis equal;
axis off;
box on;
%hold on
%annotation(gca,'line',[25/cols 530/rows],[420/cols 70/rows]);

%Now set the alpha map for the NaN region
z = ifg(:,:);
z(~isnan(ifg)) = 1;
z(isnan(ifg)) = 0;
alpha(z);
set(gca, 'color', [0 0 0]);

%plot profiles
subplot(1,2,2);
prof(:,1)=prof(:,1)-90;
sumprof(:,1)=sumprof(:,1)-90;
plot(prof(:,1),prof(:,3),'.');
hold on
plot(sumprof(:,1),sumprof(:,2),'*','color','r');
axis([-100,140,-8,8]);
axis square;
grid on
hold off
