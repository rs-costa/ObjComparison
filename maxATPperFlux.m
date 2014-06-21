function FBAsolution=maxATPperFlux(glob,model,netcode)
% Maximization of ATP yield per flux unit
[nMets,nRxns]=size(model.S);
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

initArgs{1}='max';
initArgs{2}=.5;
nOpt=500;

allObjValues=zeros(nOpt,1);
allSolutions=zeros(nRxns,nOpt);
allStats=zeros(nOpt,1);
opt=optimset('Algorithm','interior-point','Display','none');

parfor i=1:nOpt
    x0=randomObjFBASol(model,initArgs);
    if isempty(x0)
        continue;
    end
    delta0plus=x0;delta0plus(delta0plus>0)=0;
    delta0minus=x0;delta0minus(delta0minus>0)=0;
    x0=[x0;delta0plus;delta0minus];
    
    [x,f,stat]=fmincon(obj,x0,A,b,Aeq,beq,lb,ub,[],opt);
    allStats(i)=stat;
    allObjValues(i)=-f;
    allSolutions(:,i)=x(1:nRxns);

end

[FBAsolution.f, idx]=max(allObjValues);
FBAsolution.x=allSolutions(:,idx);
FBAsolution.exitflag=allStats(idx);

if FBAsolution.exitflag<0
    FBAsolution.minE=-1;
else
    FBAsolution.minE=eclDistance(glob,FBAsolution.x)/100;
end
FBAsolution.maxE=FBAsolution.minE;

end
