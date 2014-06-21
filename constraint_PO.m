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