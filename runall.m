%initCobraToolbox();
%global exp_flux var_exp indices ite bmrate multiplier nRxns;
clear;
clc;
s=matlabpool('size');
if s==0
	matlabpool open 2;
end
pctRunOnAll warning('off','all');
% for task=3:4
% 	simulation(1,0,task);	
%     simulation(2,0,task);
% end
for task=[1,2,4]
    simulation(1,0,task);
    simulation(2,0,task);
    simulation(3,0,task);
end
% simulation(1,0,4);	
% simulation(2,0,4);
matlabpool close;
