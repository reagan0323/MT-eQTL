# Multi-Tissue eQTL 

We provide Matlab codes for multi-tissue cis-eQTL analysis. The methods use transformed gene-SNP correlations across tissues as input, and jointly model them as a mixture Gaussian distribution. We use empirical Bayes approach to provide estimation of model parameters that capture global eQTL configurations and patterns in multiple tissues. We also provide a local false discovery rate approach to control the FDR of eQTL discoveries with different configurations. More details of the method can be found in the manuscripts https://arxiv.org/abs/1311.2948 (MT-eQTL) and https://arxiv.org/abs/1701.05426 (HT-eQTL).  

The MT-eQTL folder contains Matlab codes for model fitting and inference under the MT-eQTL model. The m file "Sim analysis.m" is a demo of a comprehensive analysis. It replicates the simulation study in the MT-eQTL manuscript at https://arxiv.org/abs/1311.2948.

The HT-eQTL folder contains Matlab codes for the HT-eQTL method, which is an extension of the MT-eQTL method to a large number of tissues. 
