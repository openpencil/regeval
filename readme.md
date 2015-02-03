#### Scope
The `regeval` repository provides R scripts that implement a simulation design for comparing a suite of regression methods for high-dimensional microbiome data. For the complete background, simulation and model specifications as well as evaluation results, please review:

> Shankar J, Szpakowski S et al. A systematic evaluation of high-dimensional, ensemble-based regression for exploring large model spaces in microbiome analyses. BMC Bioinformatics. 2015 Feb;16(1):31+. Available from: http://dx.doi.org/10.1186/s12859-015-0467-6.


#### R code

The files `regeval_simulation.R` and `regeval_analysis.R` are the entry points for the running the scripts.

* `regeval_simulation.R` implements the simulation design and provides evaluation and graphing routines for a systematic comparison of approaches.
* `regeval_analysis.R` applies the regression approaches on the example data and provides graphing routines for comparing the findings from the approaches.

Please step-through the code and comments within the following R scripts for detailed instructions.

| File       | Description  | Type |
|:----------|:-----------------|:----|
|`example_dataset.rda`    | An example design matrix| Data |
|`example_response.rda` | An example response vector| Data |
|`regeval_packages.R`| Installs all packages needed for the evaluation and loads the libraries | Libraries |
|`regeval_algorithms.R` | All the regression algorithms used in the evaluation. | Functions |
|`regeval_simulation.R`| Implementation of the simulation design and evaluation for a systematic comparison| Simulation |
|`regeval_graphing.R`| Graphing routines for data generated from evaluation. | Evaluation + Visualization |
|`regeval_analysis.R` | Application of the algorithms on the example data + Comparison of findings | Analysis + Visualization|
|`regeval_colorlegend.R` | Corrplot color legend | Visualization |
|`regeval_corrplot.R`| Slightly modified corrplot code |  Visualization |
|`regeval_colored_dendrogram.R` | Slightly modified cluster dendrogram code | Visualization |
