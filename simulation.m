
function simulation(md,env,task)
    gr_start=1;
%     if md==0 
%         gr_end=6;
%         envname={'ishii0.1C' 'ishii0.4C' 'ishii0.7C' 'holm0.67B' 'perrenoud0.65B' 'emmerling0.09N'};
%         dilution=[0.1 0.4 0.7 0.67 0.65 0.09];
%     else       
%         gr_end=12;
%         envname={'nanchen0.09C' 'nanchen0.4C' 'emmerling0.09N' 'emmerling0.09C' 'emmerling0.4C' 'yang0.1C' 'yang0.55C' 'ishii0.1C' 'ishii0.4C' 'ishii0.7C' 'perrenoud0.65B' 'holm0.67B'};
%         dilution=[0.09 0.4 0.09 0.09 0.4 0.1 0.55 0.1 0.4 0.7 0.65 0.67];
%     end
    gr_end=12;
    envname={'nanchen0.09C' 'nanchen0.4C' 'emmerling0.09N' 'emmerling0.09C' 'emmerling0.4C' 'yang0.1C' 'yang0.55C' 'ishii0.1C' 'ishii0.4C' 'ishii0.7C' 'perrenoud0.65B' 'holm0.67B'};
    dilution=[0.09 0.4 0.09 0.09 0.4 0.1 0.55 0.1 0.4 0.7 0.65 0.67];
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
    %           3-Genome scale model iAF1260
    %%% env: experimental data, if not specified, all are taken
    %%% task:   1-default. Use glucose uptake+additional constraints.
    %       2-pairs of those constraints.
    %		3-run all linear objective function and plot fidelity points (not bars).
    %       4-scatter plot
    %       5-run multiple objective functions optimization with 2
    %       objectives: MaxBM and other linear objectives.
    %       6-run multiple objective functions optimization with 3
    %       objectives: MaxBM, MaxATP, and other linear objectives.
    %       7-run multiple objective functions optimization with 3
    %       objectives: MaxBM, Min Flux, and other linear objectives.
    %       8-run MO for pairs and tuple in MaxBM, MaxATP, minFlux
    mdname={'ec_core_model' 'ec_schuetz_model' 'ec_iaf1260_model'};
    
    taskname={'glucose+coarse_grained' 'pairsOfConstraint' 'allLinearOFs' 'scatterPlot' 'biomass+X' 'biomass+ATP+X' 'biomass+minFlux+X' 'biomass+ATP+minFlux'};    
    load_model=strcat('models/',mdname{md});  
    load_data=strcat('data/',mdname{md},'.mat');
    
    glob.ite=100;
    tit={'max BM','max ATP','min Flux','max BM/flux','max ATP/flux','min Rd','min ATPprod','max ATPprod'};

    
%     noOfObj=length(tit)
%     log=zeros(30,4*(noOfObj+1)); % +growth rate (biomass), EX_glc, EX_o2
 
    for gr=gr_start:gr_end
        savedir=strcat('conf-Copy/',mdname{md},'/',envname{gr});
        mkdir(savedir);
        bmrate=dilution(gr);
        %%% Load model and experimental data
        load(load_model);
        load(load_data);
% 	%For paper only
% 	if md~=3
% 		strainM=strainM(:,[8,9,10,12,11,3]);
% 		strainM_var=strainM_var(:,[8,9,10,12,11,3]);
% 		idx=idx(:,[8,9,10,12,11,3]);
% 	end
	%%%%%%%%%%%%%%%
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
        elseif md==3
            multiplier=strainM(1005,gr)/bmrate;
        end
    
%         log(:,1+(1+noOfObj)*(gr-1))=[glob.exp_flux;dilution(gr);-consump_rate(1);-consump_rate(2)];
        
        glob.indices=index;
    
        for i=1:length(model.mets)
            model.ub(model.ub~=10000)=model.ub(model.ub~=10000)*multiplier;
            model.lb(model.lb~=-10000)=model.lb(model.lb~=-10000)*multiplier;
        end
        if md==1 || md==3
            % No  constraint
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
        elseif task==5
            bm={'Biomass_Ecoli_core_w_GAM' 'biomass' 'Ec_biomass_iAF1260_core_59p81M'};
            % range of biomass flux
            model = changeObjective(model,bm{md});
            fba= optimizeCbModel(model,'min');                
            bmMin= fba.f;
            fba= optimizeCbModel(model,'max');                
            bmMax= fba.f;
            
            tit={};
            solution=[];
            for i=1:length(model.mets)
                if (md==1 || md==3) && ~isempty(strfind(model.mets{i},'[e]'))  
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
                
                mindist=pareto2(glob,amodel,@maxX,md,task,bmMin,bmMax);              
                tit=[tit strcat('+Max',{' '},model.mets{i})];
                solution=[solution;mindist];
%                 fprintf('model: %s, dataset: %s, metabolite: %d/%d\n',mdname{md},envname{gr},i,length(model.mets));
            end
            
            mindist=pareto2(glob,model,@minRedox,md,task,bmMin,bmMax);
            solution=[solution;mindist];
            tit=[tit '+Min Rd'];

            mindist=pareto2(glob,model,@minATPProd,md,task,bmMin,bmMax);    
            solution=[solution;mindist];
            tit=[tit '+Min ATPprod'];

            mindist=pareto2(glob,model,@maxATPProd,md,task,bmMin,bmMax);    
            solution=[solution;mindist];
            tit=[tit '+Max ATPprod'];

            mindist=pareto2(glob,model,@minManhattan,md,task,bmMin,bmMax);    
            solution=[solution;mindist];
            tit=[tit '+Min Flux']; 
            
            tit=tit(solution>0);
            solution=solution(solution>0);

            [solution,id]=sort(solution);
            tit=tit(id);
	    
        elseif task==6
            bm={'Biomass_Ecoli_core_w_GAM' 'biomass' 'Ec_biomass_iAF1260_core_59p81M'};
            atpm={'ATPM' 'maint' 'ATPM'};
            % range of biomass flux
            model = changeObjective(model,bm{md});
            fba= optimizeCbModel(model,'min');                
            bmMin= fba.f;
            fba= optimizeCbModel(model,'max');                
            bmMax= fba.f;
            % range of ATPM 
            model = changeObjective(model,atpm{md});
            fba= optimizeCbModel(model,'min');                
            atpMin= fba.f;
            fba= optimizeCbModel(model,'max');                
            atpMax= fba.f;
            tit={};
            solution=[];
            for i=1:length(model.mets)
                if (md==1 || md==3) && ~isempty(strfind(model.mets{i},'[e]'))  
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
                
                mindist=pareto3(glob,amodel,@maxX,md,task,bmMin,bmMax,atpMin,atpMax);              
                tit=[tit strcat('+Max',{' '},model.mets{i})];
                solution=[solution;mindist];
%                 fprintf('model: %s, dataset: %s, metabolite: %d/%d\n',mdname{md},envname{gr},i,length(model.mets));
            end
            
            mindist=pareto3(glob,model,@minRedox,md,task,bmMin,bmMax,atpMin,atpMax);
            solution=[solution;mindist];
            tit=[tit '+Min Rd'];

            mindist=pareto3(glob,model,@minATPProd,md,task,bmMin,bmMax,atpMin,atpMax);    
            solution=[solution;mindist];
            tit=[tit '+Min ATPprod'];

            mindist=pareto3(glob,model,@maxATPProd,md,task,bmMin,bmMax,atpMin,atpMax);    
            solution=[solution;mindist];
            tit=[tit '+Max ATPprod'];

            mindist=pareto3(glob,model,@minManhattan,md,task,bmMin,bmMax,atpMin,atpMax);    
            solution=[solution;mindist];
            tit=[tit '+Min Flux']; 
            
            tit=tit(solution>0);
            solution=solution(solution>0);

            [solution,id]=sort(solution);
            tit=tit(id);
	    
              
        elseif task==7
            bm={'Biomass_Ecoli_core_w_GAM' 'biomass' 'Ec_biomass_iAF1260_core_59p81M'};
            % range of biomass flux
            model = changeObjective(model,bm{md});
            fba= optimizeCbModel(model,'min');                
            bmMin= fba.f;
            fba= optimizeCbModel(model,'max');                
            bmMax= fba.f;

            tit={};
            solution=[];
            for i=1:length(model.mets)
                if (md==1 || md==3) && ~isempty(strfind(model.mets{i},'[e]'))  
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
                max_sol=optimizeCbModel(model,'max');
                min_sol=optimizeCbModel(model,'min');
                min_X=0;
                max_X=10000;
                if max_sol.stat~=1
                    max_X=max_sol.f;
                end
                if min_sol.stat~=1
                    min_X=min_sol.f;
                end
                                
                mindist=pareto3_2(glob,amodel,i,md,task,bmMin,bmMax,min_X,max_X);              
                tit=[tit strcat('+Max',{' '},model.mets{i})];
                solution=[solution;mindist];
%                 fprintf('model: %s, dataset: %s, metabolite: %d/%d\n',mdname{md},envname{gr},i,length(model.mets));
            end
            
            tit=tit(solution>0);
            solution=solution(solution>0);

            [solution,id]=sort(solution);
            tit=tit(id);
	     elseif task==8
            bm={'Biomass_Ecoli_core_w_GAM' 'biomass' 'Ec_biomass_iAF1260_core_59p81M'};
            atpm={'ATPM' 'maint' 'ATPM'};
            % range of biomass flux
            model = changeObjective(model,bm{md});
            fba= optimizeCbModel(model,'min');                
            bmMin= fba.f;
            fba= optimizeCbModel(model,'max');                
            bmMax= fba.f;
            % range of ATPM 
            model = changeObjective(model,atpm{md});
            fba= optimizeCbModel(model,'min');                
            atpMin= fba.f;
            fba= optimizeCbModel(model,'max');                
            atpMax= fba.f;
            tit={};
            solution=[];
            
            %MaxBM+MaxATP           
            mindist=pareto2(glob,model,@maxATP,md,task,bmMin,bmMax);    
            solution=[solution;mindist];
            tit=[tit 'MaxBM+MaxATP'];
            %MaxBM+MinFLux           
            mindist=pareto2(glob,model,@minManhattan,md,task,bmMin,bmMax);    
            solution=[solution;mindist];
            tit=[tit 'MaxBM+MinFLux']; 
            %MaxATP+MinFLux           
            mindist=pareto2(glob,model,@minManhattan,md,task,atpMin,atpMax);    
            solution=[solution;mindist];
            tit=[tit 'MaxATP+MinFLux']; 
            %MaxBM+MaxATP+MinFlux
            mindist=pareto3(glob,model,@minManhattan,md,task,bmMin,bmMax,atpMin,atpMax);    
            solution=[solution;mindist];
            tit=[tit 'MaxBM+MaxATP+MinFlux']; 


            
            tit=tit(solution>0);
            solution=solution(solution>0);

            [solution,id]=sort(solution);
            tit=tit(id);
	          
        end
    
        fname=strcat(savedir,'/',taskname{task});
        generateResults(glob,solution,tit,constit,fname,task);
        fprintf('----------FINISH TASK %d:  %s, %s----------\n',task,mdname{md},envname{gr});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         for k=1:noOfObj
%             log(:,1+(noOfObj+1)*(gr-1)+k)=solution{1}{1,k}.xmin;
%         end


    end
%    matlabpool close;
%     dlmwrite(strcat('results',mdname{md},'\','tmp.log'),log);
 end
