function FBAsolution=maxATP(glob,model,netcode,taskcode,minNorm)
    if exist('minNorm','var')
        if isempty(minNorm)
            minNorm=true;
        end
    else
        minNorm=true;
    end
    [nMets,nRxns] = size(model.S);
    model.c=zeros(nRxns,1);
    if netcode==1 || netcode ==3
        model=changeObjective(model,'ATPM');
    elseif netcode==2
        model=changeObjective(model,'maint');
    end

    A=[model.S(model.csense=='L',:);-model.S(model.csense=='G',:)];
    b=[model.b(model.csense=='L',:);-model.b(model.csense=='G',:)];
    Aeq=model.S(model.csense=='E',:);
    beq=model.b(model.csense=='E',:);

    if minNorm
        FBAsolution=optimizeCbModel(model,'max','one');
        if FBAsolution.stat < 1
            warning('No solution found!');
            FBAsolution.minE=-1;FBAsolution.maxE=-1;
            return;
        else
            FBAsolution.minE=eclDistance(glob,FBAsolution.x)/100;
        end
        FBAsolution.maxE=FBAsolution.minE;
    else
        FBAsolution=optimizeCbModel(model,'max');
        if FBAsolution.stat < 1
            warning('No solution found!');
            return;
        end
        x0=FBAsolution.x;
        f0=FBAsolution.f;

        Aeq=[Aeq;model.c'];
        beq=[beq;f0];
        [FBAsolution.minE,FBAsolution.maxE,FBAsolution.x]=fidelity4linear(glob,taskcode,x0,Aeq,beq,A,b,model.lb,model.ub);
    end

end
