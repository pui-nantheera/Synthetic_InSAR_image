function [] = main_batch()

%K=0.5;
%deltafault=0;
%for depth=5:1:25
% parfile = strcat(num2str(depth),'km_','K',num2str(K),'delta',num2str(deltafault),'.conf');
% main(parfile);
%end

%depth=10;
%deltafault=0;
%for K=0.1:0.1:0.9
% parfile = strcat(num2str(depth),'km_','K',num2str(K),'delta',num2str(deltafault),'.conf');
% main(parfile);
%end

%depth=10;
%K=0.5;
%for deltafault=-18:2:18
%  parfile = strcat(num2str(depth),'km_','K',num2str(K),'delta',num2str(deltafault),'.conf');
%  main(parfile);
%end

 depth=10;
 for deltafault=-8:2:-2
   for K=0.1:0.1:0.9
     parfile = strcat(num2str(depth),'km_','K',num2str(K),'delta',num2str(deltafault),'.conf');
     main(parfile);
   end
 end
