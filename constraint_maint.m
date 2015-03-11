% Exploring the cellular objective in flux balance constraint-based models.

% Author: Rafael Costa
% Affiliation: Instituto de Engenharia de Sistemas e Computadores, Investigacão e Desenvolvimento (INESC-ID), Lisboa
% and
% Center for Intelligent Systems, LAETA, IDMEC, IST, University of Lisbon 
% Author: Nguyen Hoang Son

function maintmodel=constraint_maint(model,netcode,multiplier) 
    if netcode==1 || netcode==3
        maintmodel=changeRxnBounds(model,'ATPM',8.39*multiplier,'l');
    elseif netcode==2
        maintmodel=changeRxnBounds(model,'maint',8.39*multiplier,'l');
    end
end