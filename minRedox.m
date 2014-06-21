function FBAsolution=minRedox(glob,model,netcode,taskcode,minNorm)
%Min redox potential
    if exist('minNorm','var')
        if isempty(minNorm)
            minNorm=true;
        end
    else
        minNorm=true;
    end
    if netcode==1
        [modelIrrev, matchRev, rev2irrev, irrev2rev]= convertToIrreversible(model);
        nadhId=findMetIDs(modelIrrev,'nadh[c]');
        nadphId=findMetIDs(modelIrrev,'nadph[c]');
        bmId=findRxnIDs(modelIrrev,'Biomass_Ecoli_core_w_GAM');

    %     obj_rnxs={'GAPD_f','PDH','ME2','AKGDH','MDH_f','NADTRHD','LDH_D_f','ACALD_f','ALCD2x_f','ME2','G6PDH2r_f','GND','ICDHyr_f','THD2','FRD7','SUCDi'};  
    elseif netcode==2
        [modelIrrev, matchRev, rev2irrev, irrev2rev]= convertToIrreversible(model);
        nadhId=findMetIDs(modelIrrev,'NADH');
        nadphId=findMetIDs(modelIrrev,'NADPH');
        bmId=findRxnIDs(modelIrrev,'biomass');    
    %     obj_rnxs={'gapA_f','aceEF','maeA','sucAB','mdh_f','udhA','fdhF','fdoGHI','fdnGHI_r2','ldhA_f'...
    %                'adhE_b','mhpF_b','adhP_b','adhC_b','maeB','zwf_f','gnd','icd_f','pntAB','frdABCD_f','sdhAB_r','dld','sdhABCD_b'};
    elseif netcode==3
        [modelIrrev, matchRev, rev2irrev, irrev2rev]= convertToIrreversible(model);
        nadhId=findMetIDs(modelIrrev,'nadh[c]');
        nadphId=findMetIDs(modelIrrev,'nadph[c]');
        bmId=findRxnIDs(modelIrrev,'Ec_biomass_iAF1260_core_59p81M');    
    %     obj_rnxs={'GAPD_f','PDH','ME2','AKGDH','MDH_f','NADTRHD','LDH_D_f','ACALD_f','ALCD2x_f','ME2','G6PDH2r_f','GND','ICDHyr_f','THD2pp','FRD2','FRD3','SUCDi'};  
    end
    %%Get list of reactions which produce NADH or NADPH
    nadh_rnxs=full(modelIrrev.S);
    nadh_rnxs=nadh_rnxs(nadhId,:);
    nadh_rnxs(nadh_rnxs<0)=0;

    nadph_rnxs=full(modelIrrev.S);
    nadph_rnxs=nadph_rnxs(nadphId,:);
    nadph_rnxs(nadph_rnxs<0)=0;

    obj_rnxs=nadh_rnxs+nadph_rnxs;
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
        FBAsolution.minE=-1;
        FBAsolution.maxE=-1;
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
