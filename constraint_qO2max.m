% Exploring the cellular objective in flux balance constraint-based models.

% Author: Rafael Costa
% Affiliation: Instituto de Engenharia de Sistemas e Computadores, Investigacão e Desenvolvimento (INESC-ID), Lisboa
% and
% Center for Intelligent Systems, LAETA, IDMEC, IST, University of Lisbon 
% Author: Nguyen Hoang Son

function qO2model=constraint_qO2max(model,netcode,multiplier)
    %% qO2max = 15   
    if netcode==1 || netcode==3
        qO2model=changeRxnBounds(model,'EX_o2(e)',-15*multiplier,'l');
    elseif netcode==2
        qO2model=changeRxnBounds(model,'o2',15*multiplier,'u');
    end
end