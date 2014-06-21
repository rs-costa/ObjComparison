function FBAsolution=maxBMperFlux2(glob,model,netcode)
% Maximization of biomass yield per flux unit
[nMets,nRxns]=size(model.S);
% x0=model.lb+(model.ub-model.lb).*rand(length(model.rxns),1)
if netcode==1
    obj=@(x)(-x(13)/sum(x.^2));
    nOpt=100;
elseif netcode==2
    obj=@(x)(-x(98)/sum(x.^2));
    nOpt=100;
elseif netcode==3
    obj=@(x)(-x(1005)/sum(x.^2));
    nOpt=1000;
end

A=[model.S(model.csense=='L',:);-model.S(model.csense=='G',:)];
b=[model.b(model.csense=='L',:);-model.b(model.csense=='G',:)];
Aeq=model.S(model.csense=='E',:);
beq=model.b(model.csense=='E',:);

initArgs{1}='max';
initArgs{2}=.5;

xx=generateRand(model,nOpt,netcode);

allObjValues=zeros(nOpt,1);
allSolutions=zeros(nRxns,nOpt);
allStats=zeros(nOpt,1);
opt=optimset('Algorithm','interior-point','Display','none');

parfor i=1:nOpt
    %x0=randomObjFBASol(model,initArgs);
    %if isempty(x0)
    %    continue;
    %end
    x0=xx(:,i);
    [x,f,stat]=fmincon(obj,x0,A,b,Aeq,beq,model.lb,model.ub,[],opt);
    [x,f,stat]=fmincon(obj,x,A,b,Aeq,beq,model.lb,model.ub,[],opt);
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
