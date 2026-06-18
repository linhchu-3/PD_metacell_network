Multiscale Transcriptomic Network Models of Parkinson’s Disease at the Single Cell Level

Overview: 

This repository implements the metacell-based multiscale network analysis framemwok (MCMNA) - a framework that includes metaclell construction and characterization, co-expression network analysis, causal network analysism, and key driver analysis. 

The MCMNA framwork to a single-nucleus RNA sequencing (snRNA-seq) dataset derived from from substantia nigra (SN) samples of 23 idiopathic PD cases and 9 controls, with an average age of 81, as described in a previous publication. The study profiled a total of 315,867 high quality nuclei and identified 12 distinct cell clusters.

<img width="2550" height="3300" alt="2025_04_23_PD Network workflow_v4" src="https://github.com/user-attachments/assets/60a95f02-98ae-49fd-bf4a-f94eab2a4d40" />

Dataset availability: 
The raw data of the our snRNA-seq are available from the Gene Expression Omnibus (GEO) under the accession code GSE184950. The raw data from Lee et al. are available from GEO under the accession code GSE148434. The raw data from Martirosyan et al. dataset are available from GEO under the accession code GSE243639. The raw data from Ma et al. are available from GEO under the accession code GSE253462. 

Workflow Summary:

Name: MetaCell-based Multiscale Network Analysis (MCMNA)
Type: Multiscale network analysis workflow 
Language: R

Description: 
MCMNA overcomes the inherent sparsity of single cell data by leveraging metacell-based aggregation, enabling construction of robust gene regulatory networks. Gene co-expression network analysis identifies cell-type specific gene modules associated with disease. Integration of meta-cell based differential gene expression and a Bayesian causal network systematically identifies putative driver genes and the top drivers

Workflow step: 
1. Construct and characterize metacell from single cell
2. Generate Biologically Meaningful Correlation (BMC) Matrix
3. Generate metacell-based Multiscale Embedding Gene Co-Expression Network Analysis (MEGENA)
4. Generate filtered matrix  (prunning method)
5. Generate metacell-based Bayesian Netowkr
6. Identify key drivers
