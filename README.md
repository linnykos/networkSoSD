# Purpose

This repository contains all the data, functions, scripts to run simulations and analysis, and scripts to generate plots for the paper
"Bias-adjusted spectral clustering in multi-layer stochastic block models" by Jing Lei and Kevin Z. Lin.

# Installation

This package can be installed through `devtools` in R.

```{r}
library("devtools")
devtools::install_github("linnykos/networkSoSD", subdir = "networkSoSD")
```
The package itself depends on several packages. These include the R packages `clue`, `devtools`, `irlba`, `kernlab`, `MASS`, `Matrix`, `RSpectra` and `Seurat`, as well as the Bioconductor packages `clusterProfiler`, `org.Mmu.eg.db`. The simulations require a custom-package called `customSimulator` at https://github.com/linnykos/customSimulator, which has additional package requirements of `future`, `future.apply`, `pbapply` and `progressr`.

# Data 

The raw data used in this article can be found at https://github.com/emmetaobrien/nhpatlas/tree/v01.00 (the file lmd_expression_matrix_2014-03-06.zip, which redirects you to http://blueprintnhpatlas.org/static/download -- see “Microdissection Microarray”, where our data originates from the “all normalized LMD microarray expression” file). This is the preprocessed microarray data. We then compute the Pearson correlation of these microarray samples among the genes expressed in the brain following the steps in Liu, F., Choi, D., Xie, L., & Roeder, K. (2018). Global spectral clustering in dynamic networks. Proceedings of the National Academy of Sciences, 115(5), 927-932. Importantly, the list of genes we compute the correlation of is provided on our Github, at https://github.com/linnykos/networkSoSD/blob/master/data-raw/All_human_genes.txt

Researchers interested in the preprocessed correlation matrices themselves directly can contact Kevin Z. Lin.

# Reproducing the results

## Running the analysis

To reproduce the data analysis of the monkey rhesus data, we use the following files under the `main` folder.

* `analysis.R` runs the main analysis in Section 6. It loads in the aforementioned preprocessed RData `pnas.RData` (which contains all 10 correlation matrices, available on request) as well as the `All_human_genes.txt` file (available at https://github.com/linnykos/networkSoSD/blob/master/data-raw/All_human_genes.txt). Importantly, this script screens the genes and converts the correlation matrices into adjacency matrices as described in the paper, and then runs our debiased sum of squared spectral-clustering method on all 10 adjacency matrices to cluster the genes into 8 clusters.

* `postprocess.R` loads in the results from `analysis.R` and generates the results shown in the paper. Importantly, this script generates all the adjacency matrix plots (Figures 5, 11, 12, 13, 14), the edge-density plot in Figure 15, and Gene Ontology results in Table 1.

* `graph_plotting.R` loads in the results from `analysis.R` and generates the figures shown in Figure 1.

* `analysis_sensitivity.R` runs a similar analysis to `analysis.R`, except with different number of communities and threshold for constructing the adjacency matrix. This is the stability analysis shown in Appendix E. We run this on a server via a shell script shown in `analysis_sensitivity.sh`.

* `analysis_sensitivity_postprocess.R` plots the results from `analysis_sensitivity.R`, resulting in the plots in Figure 16.

## Running the simulations

To reproduce the simulations, we use the following files under the `simulation` folder. All of the plotting scripts use functions defined in `simulation_key.R`.

* Figure 2: `simulation_rho_simple.R` generates the simulation while `simulation_rho_simple_postprocess.R` creates the plot.

* Figure 3 and 10: `simulation_rho.R` generates the simulation while `simulation_rho_postprocess.R` creates the plot.

* Figure 4: `simulation_rho_bias.R` generates the data and creates the plot.

* Figure 6: `simulation_rho_degree.R` creates the plot, using the results from `simulation_rho.R`.

* Figure 7: `simulation_rho_3d.R` generates the data and creates the plot.

* Figure 8: `simulation_rho_switching.R` generates the simulation while `simulation_rho_switching_postprocess.R` creates the plot.

* Figure 9: `simulation_rho_switching2.R` generates the simulation while `simulation_rho_switching2_postprocess.R` creates the plot.
