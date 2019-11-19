# taskFC
Code used in Greene et al., 2019

Included code:
1. betaCalc_intraSubPPI.m performs intrasubject PPI, given task regressors and region x time matrix of activity (in Greene et al., used on 7 HCP tasks in 703 subjects, S1200 release)
2. betaCalc_interSubPPI.m performs intersubject PPI
3. ridge_wrapper_allTasks.m processes PPI betas to prepare them for use in phenotype prediction, and passes those matrices to ridgeCPM_vector.m, which performs ridge CPM (as described in Gao et al., 2019, Neuroimage, https://www.ncbi.nlm.nih.gov/pubmed/31336188; code for ridge CPM can be found at https://github.com/YaleMRRC/CPM)
4. ridge_wrapper_allTasks_permTest.m - same as (3) but with permuted phenotype measures for significance testing
5. ridge_visualization.m loads prediction outputs and identifies features selected x% of the time, calculating feature- and term-level contributions to predictive models (x = 75 in Greene et al., 2019)
6. ridgeCanonNet_high_vs_lowSynch.m splits rCPM predictive features by their inter-subject synchrony and organizes these features into 10 canonical networks, yielding a network x network matrix of contributions
7. findPPIbetaSign.m calculates mean PPI betas for each term; used in Greene et al., 2019 to separate positive and negative activation and FC
8. visualize_intersubPPI_share.m processes inter-subject PPI outputs
9. synch_vs_predictiveness_modeling.m performs regression analyses relating predictive utility to inter-subject synchrony
10. kfold_family.m splits participants into folds based on family ID, keeping families together; takes as input n x 2 matrix with col1 = subject ID, col2 = family ID
