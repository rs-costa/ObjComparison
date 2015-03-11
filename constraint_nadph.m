% Exploring the cellular objective in flux balance constraint-based models.

% Author: Rafael Costa
% Affiliation: Instituto de Engenharia de Sistemas e Computadores, Investigacão e Desenvolvimento (INESC-ID), Lisboa
% and
% Center for Intelligent Systems, LAETA, IDMEC, IST, University of Lisbon 
% Author: Nguyen Hoang Son

function nadphmodel=constraint_nadph(model,netcode)
    %% NADPH constraint
    nadphmodel=model;
    nadphoverprod=15.7*0.35; %15.7 NADPH per unit of biomass, 35% overproduction (Fig 7 Nanchen et al 2006)
    if netcode==1
            A_addd=zeros(1,size(nadphmodel.S,2));
            A_addd(68)=-(1/nadphoverprod);
            A_addd(13)=1;
            nadphmodel.S=[nadphmodel.S;A_addd];
            nadphmodel.b=[nadphmodel.b;0];
            nadphmodel.csense=[nadphmodel.csense;'E'];
            nadphmodel.mets=[nadphmodel.mets;'ps_nadh'];
            nadphmodel.metNames=[nadphmodel.metNames;'ps_nadh'];
            nadphmodel.ub(92)=0;  %TDH2 = 0
    elseif netcode==2
            A_addd=zeros(1,size(nadphmodel.S,2));
            A_addd(54)=-(1/nadphoverprod);
            A_addd(98)=1;
            nadphmodel.S=[nadphmodel.S;A_addd];
            nadphmodel.b=[nadphmodel.b;0];
            nadphmodel.csense=[nadphmodel.csense;'E'];
            nadphmodel.mets=[nadphmodel.mets;'ps_nadh'];
            nadphmodel.metNames=[nadphmodel.metNames;'ps_nadh'];
            nadphmodel.ub(53)=0;  %PntAB = 0

    elseif netcode==3
            A_addd=zeros(1,size(nadphmodel.S,2));
            A_addd(1725)=-(1/nadphoverprod);
            A_addd(1005)=1;
            nadphmodel.S=[nadphmodel.S;A_addd];
            nadphmodel.b=[nadphmodel.b;0];
            nadphmodel.csense=[nadphmodel.csense;'E'];
            nadphmodel.mets=[nadphmodel.mets;'ps_nadh'];
            nadphmodel.metNames=[nadphmodel.metNames;'ps_nadh'];
            nadphmodel.ub(2222)=0;  %TDH2pp = 0
    end
end