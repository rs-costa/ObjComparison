% Exploring the cellular objective in flux balance constraint-based models.

% Author: Rafael Costa
% Affiliation: Instituto de Engenharia de Sistemas e Computadores, Investigacão e Desenvolvimento (INESC-ID), Lisboa
% and
% Center for Intelligent Systems, LAETA, IDMEC, IST, University of Lisbon 
% Author: Nguyen Hoang Son

function boundmodel=constraint_bound(model,netcode)
    %% Bound constraint
    boundmodel=model;
    if netcode==1
        cons=[1:10,12,15,18,40:95];
    elseif netcode==2
        cons=[1:52,60:83];
    elseif netcode==3
        cons=[];
         for i=1:length(model.mets)
                if i~= 375 && i~= 1005 && isempty(strfind(model.mets{i},'[e]')) && isempty(strfind(model.mets{i},'ex')) && isempty(strfind(model.mets{i},'pp'))
                    cons=[cons i];
                end
         end
    end
    for k=1:length(cons)
            if boundmodel.ub(cons(k))==10000
                boundmodel.ub(cons(k))=200;
            end
            if boundmodel.lb(cons(k))==-10000
                boundmodel.lb(cons(k))=-200;
            end
    end
end