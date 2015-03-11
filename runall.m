% Exploring the cellular objective in flux balance constraint-based models.

% Author: Rafael Costa
% Affiliation: Instituto de Engenharia de Sistemas e Computadores, Investigacão e Desenvolvimento (INESC-ID), Lisboa
% and
% Center for Intelligent Systems, LAETA, IDMEC, IST, University of Lisbon 
% Author: Nguyen Hoang Son

%initCobraToolbox();
%global exp_flux var_exp indices ite bmrate multiplier nRxns;
clear;
clc;
s=matlabpool('size');
if s==0
	matlabpool open 2; % 2 cpus used in parallel
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
