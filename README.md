FBAcomparison
=============
FBAs with several objective functions for different E.coli models and datasets.

I. PREQUISITES

-MATLAB 2012b with Global Optimization Toolbox and Parallel Computing.

-Cobra toolbox 2.0.5 (GLPK+libSBML): http://opencobra.sourceforge.net/openCOBRA
 
II. CONTENT

1) Scripts:

-runall.m: control script to run all experiments in the paper.

-simulation.m: different FBA implementations go here.

-constraint_*.m: scripts to set up constraints for FBA.

-max/min*.m: objective functions.

2) Data:

-./data/Expdata_ec_core_model.mat: experimental datasets containing mapped fluxes values from different data sources to the Core model.

-./dataExpdata_schuetz_core_model.mat: experimental datasets containing mapped fluxes values from different data sources to the Schuetz model.

3) Models:

-./models/ec_core_model.xml: SBML format for the Core model.

-./models/ec_core_model.mat: Matlab format for the Core model.

-./models/ec_schuetz_model.xml: SBML format for the Schuetz model.

-./models/ec_schuetz_model.mat: Matlab format for the Schuetz model.

4) Additional files

-./supp/matching.xlsx: Microsoft excel file containing the matching between experimental fluxes from different data sources to reactions in 
reconstruction models.

-./supp/plots.docx: Microsoft word file containing all figures generated.

