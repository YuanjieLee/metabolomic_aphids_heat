Overview

This repository contains R scripts and datasets for analyzing metabolomics experiments descriped in manuscript entitled "Metabolite profiling reveals insight into interspecific variations in heat tolerance of three co-existing pest species". 
The scripts perform various tasks such as survival curve plotting, principal component analysis (PCA), and univariate analysis.

Contents
Scripts
1-Survive_curve_fig1.R
Script for computing median lethal time and generating survival curves to visualize experimental data.

2-exp1_pca.R
Performs PCA for Experiment 1, exploring variance and identifying patterns between different treatments.

3-exp2_pca.R
Performs PCA for Experiment 2, exploring variance and identifying patterns between different treatments.

4-Univariate_analysis-experiment1.R
Conducts univariate statistical analysis for Experiment 1 to identify changes in individual metabolites level.

5-Univariate_analysis-experiment2.R
Conducts univariate statistical analysis for Experiment 2 to identify changes in individual metabolites level.

Data

aphid_surv.csv:
Data file containing survival data for heat treatment. This is used as input for Lethal median time computating and survival curve plot.

metabo_exp1.txt:
Data file containing metabolites measurements for Experiment 1. This is used as input for the PCA analysis in 2-exp1_pca.R.

metabo_exp2.txt:
Data file containing metabolites measurements for Experiment 2. This is used as input for the PCA analysis in 3-exp2_pca.R.

metabo_average.csv and metabo_foldchange.csv:

Data file containing metabolites measurements, fold change, and averaged value by treatment for Experiment 1 and 2. This is used as input for the 4-Univariate_analysis-experiment1.R and 5-Univariate_analysis-experiment2.R
