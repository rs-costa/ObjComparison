Experimental datasets, metabolic models and MATLAB scripts from:

R S Costa, S Nguyen, A Hartman, S Vinga "Exploring the cellular objective in flux balance constraint-based models", Lecture Notes in Computer Science: 8859, 211-224 (2014). DOI:10.1007/978-3-319-12982-2_15 (URL: http://bit.ly/1DDdwdV)



ObjComparison
=============
Systematic comparison of several objective functions for different E. coli models and datasets.

I. PREQUISITES

-MATLAB 2012b with Global Optimization Toolbox and Parallel Computing.

-Cobra toolbox 2.0.5 (GLPK+libSBML): http://opencobra.sourceforge.net/openCOBRA
 
II. CONTENT

1) Scripts:

-runall.m: control script to run all experiments in the paper.

-simulation.m: different FBA implementations go here.

-constraint_*.m: scripts to set up constraints for FBA.

-max/min*.m: objective functions.

-pareto*.m: pareto optimization with 2 or 3 objective functions. The algorithm being applied is epsilon-constraint [1]. Not used in the article.

2) Data:

-./data/Expdata_ec_core_model.mat: experimental datasets containing mapped fluxes values from different data sources to the Core model.

-./dataExpdata_ec_schuetz_model.mat: experimental datasets containing mapped fluxes values from different data sources to the Schuetz model. Not in use for the article.

-./dataExpdata_ec_iaf1260_model.mat: experimental datasets containing mapped fluxes values from different data sources to the iAF1260 genome-scale model.

3) Models:

-./models/ec_core_model.mat: Matlab format for the Core model [2].

-./models/ec_schuetz_model.mat: Matlab format for the Schuetz model [3]. Not used in the article.

-./models/ec_iaf1260_model.mat: Matlab format for the Genome-scale model [4].

4) Additional files

-./supp/matching.xlsx: Microsoft excel file containing the matching between experimental fluxes from different data sources to reactions in 
reconstruction models.

-./supp/plots.docx: Microsoft word file containing all figures generated.


REFERENCES

[1] Jahn, J.: Vector Optimization: Theory, Applications, and Extensions (Springer, Heidelberg, Germany, 2004).

[2] Orth, J.D, Fleming, R.M.T., Palsson, B.O.: Reconstruction and use of microbial metabolic networks: the core Escherichia coli metabolic model as an Educational Guide. In: Escherichia coli and Salmonella: Cellular and Modelcular Biology, ASM Press, edition 2010.

[3] Schuetz, R., Kuepfer, L., Sauer, U.: Systematic evaluation of objective functions for predicting intracellular fluxes in Escherichia coli. Molecular Systems Biology 3, 119 (2007)

[4] Feist, A.M., Henry, C.S., Reed, J.L., Krunmenacker, M., Joyce, A.R., Karp, P.D. et al: A genome-scale metabolic reconstruction for Escherichia coli K-12 MG1655 that accounts for 1260 ORFs and thermodynamic information. Molecular Systems Biology, 3, 121(2007).
