%initCobraToolbox();
%global exp_flux var_exp indices ite bmrate multiplier nRxns;
clear;
clc;
s=matlabpool('size');
ncpu=8; %number of CPU for parallel running
if s==0
	matlabpool open ncpu;
end
pctRunOnAll warning('off','all');
% for task=[1,2,4,8]
% 	simulation(1,0,task);	
%     	simulation(2,0,task);
%	simulation(3,0,task);
% end
for task=[1,2,4]
    simulation(1,0,task);
    simulation(2,0,task);
    simulation(3,0,task);
end

matlabpool close;
