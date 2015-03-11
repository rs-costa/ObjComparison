% Exploring the cellular objective in flux balance constraint-based models.

% Author: Rafael Costa
% Affiliation: Instituto de Engenharia de Sistemas e Computadores, Investigacão e Desenvolvimento (INESC-ID), Lisboa
% and
% Center for Intelligent Systems, LAETA, IDMEC, IST, University of Lisbon 
% Author: Nguyen Hoang Son

function FBAsolution=maxBMperFlux2(glob,model,netcode)
% Maximization of biomass yield per flux unit
[nMets,nRxns]=size(model.S);
if netcode==1
    obj=@(x)(-x(13)/sum(x.^2));
%     md0=changeObjective(model,'Biomass_Ecoli_core_w_GAM');
%     load('data/gs1_initPoints.mat');
elseif netcode==2
    obj=@(x)(-x(98)/sum(x.^2));
%     md0=changeObjective(model,'biomass');
%     load('data/gs2_initPoints.mat');
elseif netcode==3
    obj=@(x)(-x(1005)/sum(x.^2));
%     md0=changeObjective(model,'Ec_biomass_iAF1260_core_59p81M');
%     load('data/gs3_initPoints.mat');
end

A=[model.S(model.csense=='L',:);-model.S(model.csense=='G',:)];
b=[model.b(model.csense=='L',:);-model.b(model.csense=='G',:)];
Aeq=model.S(model.csense=='E',:);
beq=model.b(model.csense=='E',:);
%%Play with MultiStart or GlobalSearch. set CustomStartPointSet
% sol0=optimizeCbModel(md0,'max','one');
% if sol0.stat ~=1
	x0=generateRand(model,1,netcode,1);
	x0=x0';
% else
% 	x0=sol0.x;
% end
gs_pts=generateRand(model,glob.ite,netcode,1);

opt=optimset('Algorithm','interior-point','Display','none');
problem=createOptimProblem('fmincon','objective',obj,'x0',x0,'Aeq',Aeq,'Aineq',A,'beq',beq,'bineq',b,'lb',model.lb,'ub',model.ub,'options',opt);
ms=MultiStart('UseParallel','always');

% gs_pts=[x0';gs_pts];
tpoints=CustomStartPointSet(gs_pts);
[sol,fval,exitflag,output]=run(ms,problem,tpoints);

FBAsolution.f=-fval;
FBAsolution.x=sol;
FBAsolution.exitflag=exitflag;
FBAsolution.output=output;

if FBAsolution.exitflag<0
    FBAsolution.minE=-1;
else
    FBAsolution.minE=eclDistance(glob,FBAsolution.x)/100;
end
FBAsolution.maxE=FBAsolution.minE;
%FBAsolution.midE=FBAsolution.minE;
end
