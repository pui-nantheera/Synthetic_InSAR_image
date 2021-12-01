function [] = pcmc_batch()

for i=0:10:50
 simsliprate_horz = i;
 parfile = strcat('pcmc_',num2str(simsliprate_horz),'.conf');
 pcmc_syn(parfile);
end

%for i=0:10:50
%  pcmc_inv(i);
%end
