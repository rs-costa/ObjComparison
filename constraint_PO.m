% Exploring the cellular objective in flux balance constraint-based models.

% Author: Rafael Costa
% Affiliation: Instituto de Engenharia de Sistemas e Computadores, Investigacão e Desenvolvimento (INESC-ID), Lisboa
% and
% Center for Intelligent Systems, LAETA, IDMEC, IST, University of Lisbon 
% Author: Nguyen Hoang Son

function POmodel=constraint_PO(model,netcode)
    %% P/O=1    
    if netcode==1
        POmodel=changeRxnBounds(model,'NADH16',0,'b');
        POmodel=changeRxnBounds(POmodel,'CYTBD',0,'b');
    elseif netcode==2
        POmodel=changeRxnBounds(model,'nuo',0,'b');
        POmodel=changeRxnBounds(POmodel,'cydAB',0,'b');
    elseif netcode==3
        POmodel=changeRxnBounds(model,'NADH16pp',0,'b');
        POmodel=changeRxnBounds(POmodel,'CYTBDpp',0,'b');
    end
end