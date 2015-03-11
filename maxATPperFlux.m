% Exploring the cellular objective in flux balance constraint-based models.

% Author: Rafael Costa
% Affiliation: Instituto de Engenharia de Sistemas e Computadores, Investigacão e Desenvolvimento (INESC-ID), Lisboa
% and
% Center for Intelligent Systems, LAETA, IDMEC, IST, University of Lisbon 
% Author: Nguyen Hoang Son

function FBAsolution=maxATPperFlux(glob,model,netcode)
% Maximization of ATP yield per flux unit
atpm={'ATPM' 'maint' 'ATPM'};
atpRxns=atpm{netcode};
model = changeObjective(model,atpRxns);
fba= optimizeCbModel(model,'min');                
atpMin= fba.f;
fba= optimizeCbModel(model,'max');                
atpMax= fba.f;

sol0=optimizeCbModel(model,'max','one');
if sol0.stat ~=1
	x0=generateRand(model,1,netcode,2);
	x0=x0';
else
	x0=sol0.x;
end
nStep=100;
nVar=size(model.S,2);
    
H=sparse(1:glob.nRxns,1:glob.nRxns,ones(glob.nRxns,1),nVar,nVar);
opt=optimset('Algorithm','interior-point-convex','Display','none');
f=zeros(nVar,1);

xValues=[]; objValues=[]; flagValues=[];
parfor i=0:nStep
    tmp=atpMin+i*(atpMax-atpMin)/nStep;
    md=changeRxnBounds(model,atpRxns,tmp,'b');
   
    A=[md.S(md.csense=='L',:);-md.S(md.csense=='G',:)];
    b=[md.b(md.csense=='L',:);-md.b(md.csense=='G',:)];
    Aeq=md.S(md.csense=='E',:);
    beq=md.b(md.csense=='E',:);
    [sol,fval,exitflag]=quadprog(H,f,A,b,Aeq,beq,md.lb,md.ub,x0,opt);
    
    if exitflag > 0
        objValues=[objValues tmp/(2*fval)];
        xValues= [xValues sol];
        flagValues= [flagValues exitflag];
    end

end
% fprintf('%.5f =============> %.5f\n',atpMin,atpMax);
[val, idx]=max(objValues);
FBAsolution.f=val;
FBAsolution.x=xValues(:,idx);
FBAsolution.exitflag=flagValues(idx);

if isempty(FBAsolution.x) || FBAsolution.exitflag<=0
    FBAsolution.minE=-1;
else
    FBAsolution.minE=eclDistance(glob,FBAsolution.x)/100;
end
FBAsolution.maxE=FBAsolution.minE;
end
