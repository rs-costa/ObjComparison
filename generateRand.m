function xx=generateRand(model,numpoint,netcode,funcode)
    if netcode==1 || netcode==3
        % No more additional constraints
        model=changeRxnBounds(model,'ATPM',0,'l');
        model=changeRxnBounds(model,'EX_glc(e)',-100,'b');
    elseif netcode ==2
        model=changeRxnBounds(model,'maint',0,'l');
        model=changeRxnBounds(model,'ptsGHI',100,'b');

    end
    [nMets,nRxns]=size(model.S);
    xx=zeros(numpoint,nRxns);
    %% Old implementation (without option funcode)
% 	initArgs{1}='max';
% 	initArgs{2}=.5;
% 	parfor i=1:numpoint
% 		xx(i,:)=randomObjFBASol(model,initArgs);	
% 	end		

    %% New implementation from here
    if funcode==1 % biomass per flux
        bm={'Biomass_Ecoli_core_w_GAM' 'biomass' 'Ec_biomass_iAF1260_core_59p81M'};
        rxns=bm{netcode};
    elseif funcode ==2 %atp per flux
        atpm={'ATPM' 'maint' 'ATPM'};
        rxns=atpm{netcode};
    end
    
    nVar=nRxns; % funcode involve...
    
    H=sparse(1:nRxns,1:nRxns,ones(nRxns,1),nVar,nVar);
    opt=optimset('Algorithm','interior-point-convex','Display','none');
    f=zeros(nVar,1);
    
    model = changeObjective(model,rxns);
    fba= optimizeCbModel(model,'min');                
    rMin= fba.f;
    fba= optimizeCbModel(model,'max');                
    rMax= fba.f;
    x0=fba.x;
    parfor i=1:numpoint
        tmp=rMin+i*(rMax-rMin)/numpoint;
        md=changeRxnBounds(model,rxns,tmp,'b');

        A=[md.S(md.csense=='L',:);-md.S(md.csense=='G',:)];
        b=[md.b(md.csense=='L',:);-md.b(md.csense=='G',:)];
        Aeq=md.S(md.csense=='E',:);
        beq=md.b(md.csense=='E',:);
        [sol,fval,exitflag]=quadprog(H,f,A,b,Aeq,beq,md.lb,md.ub,x0,opt);
        if exitflag > 0
            xx(i,:)= sol;
        end

    end

end
