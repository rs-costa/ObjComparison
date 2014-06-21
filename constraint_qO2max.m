function qO2model=constraint_qO2max(model,netcode,multiplier)
    %% qO2max = 15   
    if netcode==1 || netcode==3
        qO2model=changeRxnBounds(model,'EX_o2(e)',-15*multiplier,'l');
    elseif netcode==2
        qO2model=changeRxnBounds(model,'o2',15*multiplier,'u');
    end
end