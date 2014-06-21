function FBAsolution=maxATPPerFluxUnit(glob,model,netcode)
% Maximization of ATP yield per flux unit
[nMets,nRxns]=size(model.S);
x0=ones(3*nRxns,1);
% x0=model.lb+(model.ub-model.lb).*rand(length(model.rxns),1)
if netcode==1
    obj=@(x)(-x(11)/sum(x((nRxns+1):3*nRxns)));
elseif netcode==2
    obj=@(x)(-x(97)/sum(x((nRxns+1):3*nRxns)));
elseif netcode==3
    obj=@(x)(-x(375)/sum(x((nRxns+1):3*nRxns)));
end

lb = [model.lb;zeros(2*nRxns,1)];
ub = [model.ub;10000*ones(2*nRxns,1)];

A=[ -speye(nRxns,nRxns) -speye(nRxns,nRxns) sparse(nRxns,nRxns);
    speye(nRxns,nRxns) sparse(nRxns,nRxns) -speye(nRxns,nRxns)];
b=zeros(2*nRxns,1);
Aeq=[model.S sparse(nMets,2*nRxns)];
beq=model.b;


opt=optimset('Algorithm','interior-point','Display','none');
problem=createOptimProblem('fmincon','objective',obj,'x0',x0,'Aeq',Aeq,'Aineq',A,'beq',beq,'bineq',b,'lb',lb,'ub',ub,'options',opt);
ms=MultiStart('UseParallel','always');
[sol,fval,exitflag,output]=run(ms,problem,glob.ite);

FBAsolution.f=-fval;
FBAsolution.x=sol;
FBAsolution.exitflag=exitflag;
FBAsolution.output=output;

FBAsolution.minE=eclDistance(glob,FBAsolution.x)/100;
FBAsolution.maxE=FBAsolution.minE;


end
