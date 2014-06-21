function maintmodel=constraint_maint(model,netcode,multiplier) 
    if netcode==1 || netcode==3
        maintmodel=changeRxnBounds(model,'ATPM',8.39*multiplier,'l');
    elseif netcode==2
        maintmodel=changeRxnBounds(model,'maint',8.39*multiplier,'l');
    end
end