function mindist=pareto3_2(glob,model,met,netcode,taskcode,bmMin,bmMax,minX,maxX)
    bm={'Biomass_Ecoli_core_w_GAM' 'biomass'};
    nstep=100;
    bmEps= (bmMax-bmMin)/nstep;    
    xEps=(maxX-minX)/nstep;
%     f=zeros(1,3);    
%     ftmp=0;
    dist=[];
    
    parfor i=0:nstep
        e1=bmMin+i*bmEps;
        %model=changeRxnBounds(model,bm{netcode},e1,'l');
        for j=0:nstep
            e2=minX+j*xEps;
            md=changeRxnBounds(model,bm{netcode},e1,'l');
            md.lb(met)=e2;
            
            solution=minManhattan(glob,md,netcode,taskcode);      
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
