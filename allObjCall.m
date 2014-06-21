function solution=allObjCall(glob,model,netcode,taskcode)
%%% Objective functions for FBA %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pFBA=false;
solution{1} = maxBiomass(glob,model,netcode,taskcode);
solution{2} = maxATP(glob,model,netcode,taskcode);
solution{3} = minManhattan(glob,model,netcode,taskcode);
solution{4} = maxBMperFlux2(glob,model,netcode);
solution{5} = maxATPperFlux2(glob,model,netcode);
% solution{4} = maxBiomassPerFluxUnit(glob,model,netcode);
% solution{5} = maxATPPerFluxUnit(glob,model,netcode);
solution{6} = minRedox(glob,model,netcode,taskcode);
solution{7} = minATPProd(glob,model,netcode,taskcode);
solution{8}= maxATPProd(glob,model,netcode,taskcode);

end
