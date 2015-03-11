% Exploring the cellular objective in flux balance constraint-based models.

% Author: Rafael Costa
% Affiliation: Instituto de Engenharia de Sistemas e Computadores, Investigacão e Desenvolvimento (INESC-ID), Lisboa
% and
% Center for Intelligent Systems, LAETA, IDMEC, IST, University of Lisbon 
% Author: Nguyen Hoang Son

function mindist=pareto3(glob,model,objfunc,netcode,taskcode,bmMin,bmMax,atpMin,atpMax)
    bm={'Biomass_Ecoli_core_w_GAM' 'biomass'};
    atpm={'ATPM' 'maint'};
    nstep=1000;
    bmEps= (bmMax-bmMin)/nstep;
    atpEps= (atpMax-atpMin)/nstep;            
%     f=zeros(1,3);    
%     ftmp=0;
    dist=[];
    
    parfor i=0:nstep
        e1=bmMin+i*bmEps;
        %model=changeRxnBounds(model,bm{netcode},e1,'l');
        for j=0:nstep
            e2=atpMin+j*atpEps;
            md=changeRxnBounds(model,bm{netcode},e1,'l');
            md=changeRxnBounds(md,atpm{netcode},e2,'l');
            
            solution=objfunc(glob,md,netcode,taskcode);      
	    if solution.stat == 1 
	    	%ftmp=solution.f;
	    	%xtmp=solution.x;
            	dist = [dist;solution.minE];
%           	f(1)=abs((ftmp-range(idx(1),1))/range(idx(1),1));
%           	f(2)=abs((transpose(cc{idx(2)})*xtmp-range(idx(2),1))/range(idx(2),1));
%           	f(3)=abs((transpose(cc{idx(3)})*xtmp-range(idx(3),1))/range(idx(3),1));
            end
	end
    end
%   f=[f;[-ftmp -e1 -e2]];
%   x=[x;xtmp'];
    mindist=min(dist(dist>0));
end
