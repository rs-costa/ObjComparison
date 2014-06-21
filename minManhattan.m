function FBAsolution=minManhattan(glob,model,netcode,taskcode,minNorm)
    if exist('minNorm','var')
        if isempty(minNorm)
            minNorm=true;
        end
    else
        minNorm=true;
    end
    [nMets,nRxns] = size(model.S);
    % Minimization l1-norm of the overall intracellular flux
    LPproblem.A = [ model.S sparse(nMets,2*nRxns);
                    speye(nRxns,nRxns) speye(nRxns,nRxns) sparse(nRxns,nRxns);
                    -speye(nRxns,nRxns) sparse(nRxns,nRxns) speye(nRxns,nRxns)];
    LPproblem.c = [zeros(nRxns,1);ones(2*nRxns,1)];
    LPproblem.lb = [model.lb;zeros(2*nRxns,1)];
    LPproblem.ub = [model.ub;10000*ones(2*nRxns,1)];
    LPproblem.b= [model.b;zeros(2*nRxns,1)];
    if ~isfield(model,'csense')
        % If csense is not declared in the model, assume that all
        % constraints are equalities.
         LPproblem.csense(1:nMets) = 'E';
    else % if csense is in the model, move it to the lp problem structure
        if length(model.csense)~=nMets,
            warning('Length of csense is invalid! Defaulting to equality constraints.')
            LPproblem.csense(1:nMets) = 'E';
        else
             LPproblem.csense = columnVector(model.csense);
        end
    end
    LPproblem.csense((nMets+1):(nMets+2*nRxns)) = 'G';
    LPproblem.csense = columnVector(LPproblem.csense);
    LPproblem.osense = 1;

    solution = solveCobraLP(LPproblem);

    FBAsolution.stat = solution.stat;
    % Store results
    if (solution.stat > 0)
        FBAsolution.x = solution.full(1:3*nRxns);
        FBAsolution.f = LPproblem.c'*solution.full(1:3*nRxns); 
    else
        FBAsolution.f = 0;
        FBAsolution.x = [];
        FBAsolution.minE=-1;
        FBAsolution.maxE=-1;
        return;
    end
    if minNorm
        FBAsolution.x=FBAsolution.x(1:nRxns);
        FBAsolution.minE=eclDistance(glob,FBAsolution.x)/100;
        FBAsolution.maxE=FBAsolution.minE;
    else
        x0=FBAsolution.x;
        A=[LPproblem.A(LPproblem.csense=='L',:);-LPproblem.A(LPproblem.csense=='G',:)];
        b=[LPproblem.b(LPproblem.csense=='L',:);-LPproblem.b(LPproblem.csense=='G',:)];
        Aeq=LPproblem.A(LPproblem.csense=='E',:);
        beq=LPproblem.b(LPproblem.csense=='E',:);
        Aeq=[Aeq;LPproblem.c'];
        beq=[beq;FBAsolution.f];

        [FBAsolution.minE,FBAsolution.maxE,FBAsolution.x]=fidelity4linear(glob,taskcode,x0,Aeq,beq,A,b,LPproblem.lb,LPproblem.ub);
        FBAsolution.x=FBAsolution.x(1:nRxns);
    end
    
end
