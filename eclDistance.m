% Exploring the cellular objective in flux balance constraint-based models.

% Author: Rafael Costa
% Affiliation: Instituto de Engenharia de Sistemas e Computadores, Investigacão e Desenvolvimento (INESC-ID), Lisboa
% and
% Center for Intelligent Systems, LAETA, IDMEC, IST, University of Lisbon 
% Author: Nguyen Hoang Son

function dist=eclDistance(glob,x)
noOfRatios=length(glob.exp_flux);
W=zeros(noOfRatios,noOfRatios);
for i=1:noOfRatios
 	W(i,i)=(1/sum(1./glob.var_exp))*(1/glob.var_exp(i));
end
dist=sqrt((x(glob.indices)-glob.exp_flux)'*W*(x(glob.indices)-glob.exp_flux));
end