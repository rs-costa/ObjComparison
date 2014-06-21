function x0=generateRand(model,numpoint,netcode)
    if netcode==1 || netcode==3
        % No more additional constraints
        model=changeRxnBounds(model,'ATPM',0,'l');
        model=changeRxnBounds(model,'EX_glc(e)',-100,'b');
    elseif netcode ==2
        model=changeRxnBounds(model,'maint',0,'l');
        model=changeRxnBounds(model,'ptsGHI',100,'b');

    end
    [nMets,nRxns]=size(model.S);
    x0=zeros(nRxns,numpoint);
    for i=1:nRxns
        model.c=zeros(nRxns,1);
        model.c(i)=1;
        max_sol=optimizeCbModel(model,'max');
        min_sol=optimizeCbModel(model,'min');
        min_X=model.ub(i);
        max_X=model.lb(i);
        if max_sol.stat~=1
            max_X=max_sol.f;
        end
        if min_sol.stat~=1
            min_X=min_sol.f;
        end
        x0(i,:)=min_X+(max_X-min_X)*rand(numpoint,1);
    end
end