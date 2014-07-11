
function simulation(md,env,task)
    gr_start=1;
    gr_end=6;
    envname={'ishii0.1C' 'ishii0.4C' 'ishii0.7C' 'holm0.67B' 'perrenoud0.65B' 'emmerling0.09N'};
    dilution=[0.1 0.4 0.7 0.67 0.65 0.09];
    if ~exist('md','var')
        md=1;
    end
    
    if exist('env','var')
        if(env < gr_start || env > gr_end)
            warning('Invalid culture code! Set to default: all');
        else
            gr_start = env;
            gr_end = env;
        end
    end
    
    if ~exist('task','var')
        task=1;
    end
    %%% Parsing input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% md:  1-default. E coli core model
    %           2-Schuetz model
    %%% env: experimental data, if not specified, all are taken
    %%% task:   1-default. Use glucose uptake+additional constraints for FBAs.
    %       2-FBAs using pairwise constraints.
    %		3-FBAs with all linear objective functions.
    %       4-Scatter plot
    mdname={'ec_core_model' 'ec_schuetz_model' 'ec_iaf1260_model'};
    
    taskname={'glucose+coarse_grained' 'pairsOfConstraint' 'allLinearOFs' 'scatterPlot'};    
    load_model=strcat('models/',mdname{md});  
    load_data=strcat('data/Expdata_',mdname{md},'.mat');
    
    glob.ite=100;
    tit={'max BM','max ATP','min Flux','max BM/flux','max ATP/flux','min Rd','min ATPprod','max ATPprod'};

    
%     noOfObj=length(tit)
%     log=zeros(30,4*(noOfObj+1)); % +growth rate (biomass), EX_glc, EX_o2
 
    for gr=gr_start:gr_end
        savedir=strcat('ver3/',mdname{md},'/',envname{gr});
        mkdir(savedir);
        bmrate=dilution(gr);
        %%% Load model and experimental data
        load(load_model);
        load(load_data);
        model.c=zeros(length(model.rxns),1);  
        glob.nRxns=size(model.S,2);
        
        index=idx(:,gr);
        index=index(index~=0);
        glob.exp_flux=strainM(index,gr);
        glob.var_exp=strainM_var(index,gr);
        if md==1
            multiplier=strainM(13,gr)/bmrate;
        elseif md==2
            multiplier=strainM(98,gr)/bmrate;
        end
    
%         log(:,1+(1+noOfObj)*(gr-1))=[glob.exp_flux;dilution(gr);-consump_rate(1);-consump_rate(2)];
        
        glob.indices=index;
    
        for i=1:length(model.mets)
            model.ub(model.ub~=10000)=model.ub(model.ub~=10000)*multiplier;
            model.lb(model.lb~=-10000)=model.lb(model.lb~=-10000)*multiplier;
        end
        if md==1
            % No more additional constraints
            model=changeRxnBounds(model,'ATPM',0,'l');
            model=changeRxnBounds(model,'EX_glc(e)',-100,'b');
        elseif md ==2
            model=changeRxnBounds(model,'maint',0,'l');
            model=changeRxnBounds(model,'ptsGHI',100,'b');
            
        end
        clear('solution','model_*');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Run simulation from here....
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %initCobraToolbox;

    	constit={}; 
        %% Glucose uptake rate + coarse grained constraints   
        if task==1
            
            constit={'No additional constraint','P/O=1','qO2 max=15','Maintenance','Bounds','NADPH','All constraints'};
            
            i=1;
            %No more additional constraints 
            mdl{i} = model;
            i=i+1;
            %% P/O=1    
            mdl{i}=constraint_PO(model,md);
            i=i+1;

            % Maximum Oxygen uptake = 15
            mdl{i} = constraint_qO2max(model,md,multiplier);
            i=i+1;
            % Maintenance constraint
            mdl{i} = constraint_maint(model,md,multiplier);
            i=i+1;

            % Bound
            mdl{i}= constraint_bound(model,md);
            i=i+1;
            % NADPH
            mdl{i}= constraint_nadph(model,md);
            i=i+1;
            % All constraints
            mdl{i}=constraint_PO(model,md);
            mdl{i} = constraint_qO2max(mdl{i},md,multiplier);
            mdl{i} = constraint_maint(mdl{i},md,multiplier);
            mdl{i}= constraint_bound(mdl{i},md);
            mdl{i} = constraint_nadph(mdl{i},md);      
%             mdl{i} = constraint_growth(mdl{i},md,bmrate);
                     
            for j=1:length(mdl)
                solution{j}=allObjCall(glob,mdl{j},md,task);
            end

        elseif task==2
            
            constit={'P/O=1 & qO2 max=15', 'P/O=1 & Maintenance','P/O=1 & Bounds', 'P/O=1 & NADPH',...
                    'qO2 max=15 & Maintenance', 'qO2 max=15 & Bounds', 'qO2 max=15 & NADPH',...
                    'Maintenance & Bounds', 'Maintenance & NADPH',  'Bounds & NADPH'};
         
           
            i=1;
            %% Create model with P/O constraint alone
            model_PO=constraint_PO(model,md);
            %P/O=1 and qO2 max = 15
            mdl{i} = constraint_qO2max(model_PO,md,multiplier);
            i=i+1;
            % P/O=1 and Maintenance
            mdl{i} = constraint_maint(model_PO,md,multiplier);
            i=i+1;
            % P/O=1 and Bound
            mdl{i}=constraint_bound(model_PO,md);
            i=i+1;
            % P/O=1 and NADPH
            mdl{i}=constraint_nadph(model_PO,md);
            i=i+1;

            %% Create model with qO2 max=15 
            model_O2=constraint_qO2max(model,md,multiplier);
            % qO2 max = 15 and Maintenance
            mdl{i} = constraint_maint(model_O2,md,multiplier);
            i=i+1;
            % qO2 max = 15 and Bound
            mdl{i}=constraint_bound(model_O2,md);
            i=i+1;
            % qO2 max = 15 and NADPH
            mdl{i}=constraint_nadph(model_O2,md);
            i=i+1;

            %% Create model with maintenance constraint alone
            model_maint=constraint_maint(model,md,multiplier);
            % Maintenance and Bound
            mdl{i}=constraint_bound(model_maint,md);
            i=i+1;
            % Maintenance and NADPH
            mdl{i}=constraint_nadph(model_maint,md);
            i=i+1;

            %% Model with bound constraint alone
            model_bound=constraint_bound(model,md);
            % Bound and NADPH
            mdl{i}=constraint_nadph(model_bound,md);
            i=i+1;
           
            for j=1:length(mdl)
                solution{j}=allObjCall(glob,mdl{j},md,task);
            end

        elseif task==3
            
            model.c=zeros(length(model.rxns),1);
            tit={};
            solution=[];
    
            for i=1:length(model.mets)
                if md==1 && ~isempty(strfind(model.mets{i},'[e]'))  
                    continue;
                end
                if md==2 && i>48
                    break
                end
                if ~isempty(strfind(lower(model.mets{i}),'atp'))
                    continue;
                end
                amodel=model;
                c=amodel.S(i,:);
                amodel.csense(i)='G';
                amodel.c=full(c');
                
                fba=maxX(glob,amodel,md,task);
                tit=[tit strcat('Max',{' '},model.mets{i})];
                solution=[solution;fba.minE];
            end
            %%% Other objective here...
%             exchange=maxExchangeFluxes(glob,model,md,task);
%             solution=[solution;exchange.minE];
%             tit=[tit 'Max Exchange Fluxes'];
            
            biomass=maxBiomass(glob,model,md,task);
            solution=[solution;biomass.minE];
            tit=[tit 'Max BM'];

            maxatp= maxATP(glob,model,md,task);
            solution=[solution;maxatp.minE];
            tit=[tit 'Max ATP'];

            minredox = minRedox(glob,model,md,task);
            solution=[solution;minredox.minE];
            tit=[tit 'Min Rd'];

            atpprod = minATPProd(glob,model,md,task);
            solution=[solution;atpprod.minE];
            tit=[tit 'Min ATPprod'];

            maxaptprod= maxATPProd(glob,model,md,task);
            solution=[solution;maxaptprod.minE];
            tit=[tit 'Max ATPprod'];

            manhattan=minManhattan(glob,model,md,task);
            solution=[solution;manhattan.minE];
            tit=[tit 'Min Flux']; 

            tit=tit(solution>0);
            solution=solution(solution>0);

            [solution,id]=sort(solution);
            tit=tit(id);
        elseif task==4       
            solution=allObjCall(glob,model,md,task);
            
        end
    
        fname=strcat(savedir,'/',taskname{task});
        generateResults(glob,solution,tit,constit,fname,task);
        fprintf('----------FINISH TASK %d:  %s, %s----------\n',task,mdname{md},envname{gr});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
 end
