% Exploring the cellular objective in flux balance constraint-based models.

% Author: Rafael Costa
% Affiliation: Instituto de Engenharia de Sistemas e Computadores, Investigacão e Desenvolvimento (INESC-ID), Lisboa
% and
% Center for Intelligent Systems, LAETA, IDMEC, IST, University of Lisbon 
% Author: Nguyen Hoang Son

function [minEucl,maxEucl,x]=fidelity4linear(glob,task,x0,Aeq,beq,A,b,lb,ub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Predict fidelitiy range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global exp_flux var_exp indices ite nRxns;
noOfRatios=length(glob.exp_flux);
W=zeros(noOfRatios,noOfRatios);
for i=1:noOfRatios
 	W(i,i)=(1/sum(1./glob.var_exp))*(1/glob.var_exp(i));
end
obj=@(x)((x(glob.indices)-glob.exp_flux)'*W*(x(glob.indices)-glob.exp_flux));
% norm_factor=sqrt(glob.exp_flux'*W*glob.exp_flux);
norm_factor=100;

%% For task 1,2
x=x0;
% if task==1 || task==2
if task==0
    opt=optimset('Algorithm','interior-point','Display','none');
    problem=createOptimProblem('fmincon','objective',obj,'x0',x0,'Aeq',Aeq,'Aineq',A,'beq',beq,'bineq',b,'lb',lb,'ub',ub,'options',opt);
    ms=MultiStart('UseParallel','always');
    [sol,fval,exitflag,output]=run(ms,problem,glob.ite);

    x=sol;
    minEucl=sqrt(fval)/norm_factor; %normalized according to Glucose uptake rate (=100)

    obj=@(x)(-(x(glob.indices)-glob.exp_flux)'*W*(x(glob.indices)-glob.exp_flux));
    problem=createOptimProblem('fmincon','objective',obj,'x0',x0,'Aeq',Aeq,'Aineq',A,'beq',beq,'bineq',b,'lb',lb,'ub',ub,'options',opt);
    ms=MultiStart('UseParallel','always');
    [sol,fval,exitflag,output]=run(ms,problem,glob.ite);

    maxEucl=sqrt(-fval)/norm_factor;
%% For task 3,4,5
else 
    nVar=size(Aeq,2);
    
    H=sparse(1:glob.nRxns,1:glob.nRxns,ones(glob.nRxns,1),nVar,nVar);
    opt=optimset('Algorithm','interior-point-convex','Display','none');
    f=zeros(nVar,1);
    [sol,fval,exitflag]=quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,opt);
    %l2norm=@(x)(x'*H*x);
    %opt=optimset('Algorithm','interior-point','Display','none');
    %problem=createOptimProblem('fmincon','objective',l2norm,'x0',x0,'Aeq',Aeq,'Aineq',A,'beq',beq,'bineq',b,'lb',lb,'ub',ub,'options',opt);
    %ms=MultiStart('UseParallel','always');
    %[sol,fval,exitflag,output]=run(ms,problem,glob.ite);
    x=sol;
    if exitflag >= 0
       minEucl=sqrt(obj(x))/norm_factor;
       maxEucl=minEucl;
    else
        minEucl=-1;maxEucl=-1;x=x0;
    end
end
end
