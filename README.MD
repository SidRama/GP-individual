A personalised approach for identifying disease-relevant pathways in heterogeneous diseases
===========================================================================================
This repository contains the Matlab scripts used in our [paper](https://doi.org/10.1101/738062).

Overview
---------------
Analysing data from complex diseases in a personalised manner to identify disrupted pathways can improve elucidation of the disease progression. We have developed a personalised statistical method that robustly models gene expression time-course datasets from complex heterogeneous diseases and summarises the gene-level results on the level of pathways, i.e. groups of genes involved in particular biological processes that drive the disease. By analysing three Type 1 Diabetes datasets using this method, we demonstrate that our personalised method reveals more insight into the biological processes involved in disease progression over time and in specific time intervals, than non-personalised methods. With its robust capabilities of identifying numerous disease-relevant pathways, this method could be further developed for predicting events in the progression of heterogeneous diseases, pursuing preventive treatments and even biomarker identification.
![Overview of the method. A schematic illustration of our personalised approach and a population-wide approach (combined method). In the personalised approach, we identify DEGs independently for each case-control pair and combine results at the pathway-level.](./images/overview.png)

#### Our personalised method
We present a method that models time-course data in a personalised manner using Gaussian processes in order to identify differentially expressed genes (DEGs); and combines the DEG lists on a pathway-level using a permutation-based empirical hypothesis testing in order to overcome gene-level variability and inconsistencies prevalent to datasets from heterogenous diseases. Our method can be applied to study the time-course dynamics as well as specific time-windows of heterogeneous diseases.

Prerequisites
---------------
These scripts require the following software:
- Matlab (>= r2016a)
- [GPstuff 4.7](https://research.cs.aalto.fi/pml/software/gpstuff/)

Using our method
---------------
For personalised **time-course** analysis:
1. Run `compute_ratios.m` and pass the required parameters.
	- The `probeset_file` is a path to a `.mat` file that contains the expression values for each probe-set (or gene) of a case-control pair (i.e. case expressions and control expressions separately). It also contains the time points for the case and control expressions as well as the time of seroconversion (or disease diagnosis).
2. Run `child_mapping.m` and pass the required parametets to map the differentially expressed probe-sets to genes. Perform this step only if probe-sets are used in the analysis.
3. To compute the `adjusted gemetric mean` for each pathway, run `pathway_overlap.m`.
4. Run `pathway_rand_overlap.m` to generate the null distribution as described in our [paper](https://doi.org/10.1101/738062).
5. Compute the emperical p-values for each pathway using the results from steps 4 and 5. Also, perform multiple testing correction on the p-values using the Benjamini-Hochberg procedure.

For personalised **time-window** analysis:
1. Run `compute_KL.m` and pass the required parameters.
	- The `probeset_file` is a path to a `.mat` file that contains the expression values for each probe-set (or gene) of a case-control pair (i.e. case expressions and control expressions separately). It also contains the time points for the case and control expressions as well as the time of seroconversion (or disease diagnosis).
2. Perform steps 2, 3, 4 and 5 as in the `personalised time-course analysis`.

Cite
---------------
Please cite this work as:
```
Somani, J., Ramchandran, S., & Lähdesmäki, H. (2019). A personalised approach for identifying disease-relevant pathways in heterogeneous diseases. BioRxiv, 738062.
```

Authors
---------------
- [Siddharth Ramchandran](https://www.siddharthr.com)
- [Juhi Somani](https://fi.linkedin.com/in/somani-juhi-942399120)
- [Harri Lähdesmäki](http://users.ics.aalto.fi/harrila/)

License
---------------
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.