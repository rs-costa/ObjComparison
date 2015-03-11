% Exploring the cellular objective in flux balance constraint-based models.

% Author: Rafael Costa
% Affiliation: Instituto de Engenharia de Sistemas e Computadores, Investigacão e Desenvolvimento (INESC-ID), Lisboa
% and
% Center for Intelligent Systems, LAETA, IDMEC, IST, University of Lisbon 
% Author: Nguyen Hoang Son

function FBAsolution=minATPProd(glob,model,netcode,taskcode,minNorm)
%Min ATP production
    if exist('minNorm','var')
        if isempty(minNorm)
            minNorm=true;
        end
    else
        minNorm=true;
    end
    if netcode==1
        [modelIrrev, matchRev, rev2irrev, irrev2rev]= convertToIrreversible(model);
        atpId=findMetIDs(modelIrrev,'atp[c]');
        bmId=findRxnIDs(modelIrrev,'Biomass_Ecoli_core_w_GAM');
    %     obj_rnxs={'PGK_f','PYK','SUCOAS_f','ATPS4r_f','ACKr_b'};

    elseif netcode==2
        [modelIrrev, matchRev, rev2irrev, irrev2rev]= convertToIrreversible(model);
        atpId=findMetIDs(modelIrrev,'ATP');
        bmId=findRxnIDs(modelIrrev,'biomass');
    %     obj_rnxs={'pgk_f','pykA','pykF','sucCD_f','atp'...
    %                'ackA_f','ackB_f','tdcD_f','purT_f'};    
    elseif netcode==3 
        [modelIrrev, matchRev, rev2irrev, irrev2rev]= convertToIrreversible(model);
        atpId=findMetIDs(modelIrrev,'atp[c]');
        bmId=findRxnIDs(modelIrrev,'Ec_biomass_iAF1260_core_59p81M');
    %     obj_rnxs={'PGK_f','PYK','SUCOAS_f','ATPS4rpp_f','ACKr_b'};
    end
    %%Get list of reactions which produce ATP
    obj_rnxs=full(modelIrrev.S);
    obj_rnxs=obj_rnxs(atpId,:);
    obj_rnxs(obj_rnxs<0)=0;
    obj_rnxs(bmId)=0;

    modelIrrev.csense=model.csense;
    modelIrrev.c=obj_rnxs';
    if minNorm
        FBAsolution = optimizeCbModel(modelIrrev,'min','one');
    else
        FBAsolution = optimizeCbModel(modelIrrev,'min');
    end
    if FBAsolution.stat < 1
        FBAsolution.x=[];
        FBAsolution.f=0;
        return;
    end
    x_irrev=zeros(length(rev2irrev),1);
    for i=1:length(rev2irrev)
        x_irrev(i)=sum(FBAsolution.x(rev2irrev{i}));
    end
    FBAsolution.x=x_irrev;

    if minNorm
        FBAsolution.minE=eclDistance(glob,FBAsolution.x)/100;
        FBAsolution.maxE=FBAsolution.minE;
    else
        A=[model.S(model.csense=='L',:);-model.S(model.csense=='G',:)];
        b=[model.b(model.csense=='L',:);-model.b(model.csense=='G',:)];
        Aeq=model.S(model.csense=='E',:);
        beq=model.b(model.csense=='E',:);

        x0=FBAsolution.x;
        f0=FBAsolution.f;

        Aeq=[Aeq;model.c'];
        beq=[beq;f0];

        [FBAsolution.minE,FBAsolution.maxE,FBAsolution.x]=fidelity4linear(glob,taskcode,x0,Aeq,beq,A,b,model.lb,model.ub);
    end
end
