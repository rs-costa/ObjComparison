% Exploring the cellular objective in flux balance constraint-based models.

% Author: Rafael Costa
% Affiliation: Instituto de Engenharia de Sistemas e Computadores, Investigacão e Desenvolvimento (INESC-ID), Lisboa
% and
% Center for Intelligent Systems, LAETA, IDMEC, IST, University of Lisbon 
% Author: Nguyen Hoang Son

function growthmodel=constraint_growth(model,netcode,bmrate)
    %% set growth rate as biomass flux  
    if netcode==1
        growthmodel=changeRxnBounds(model,'Biomass_Ecoli_core_w_GAM',bmrate,'b');
    elseif netcode==2
        growthmodel=changeRxnBounds(model,'biomass',bmrate,'b');
    elseif netcode==3
        growthmodel=changeRxnBounds(model,'Ec_biomass_iAF1260_core_59p81M',bmrate,'b');
    end
end