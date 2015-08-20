The Main.m file is the main file for carrying out MT-eQTL.  We assume users have two data sets: cormat.mat and df.mat.
The matrix "cormat" is of size N*K, where each entry is the Pearson correlation of a cis gene-SNP pair in one tissue. 
Each row corresponds to a cis gene-SNP pair (N pairs in total), and each column corresponds to a tissue (K tissues in total). The donors for different tissues may vary. 
The date set "df" is a K*1 vector recording degrees of freedom in all K tissues (df=the number of samples minus the number of covariates minus 2).

Running the main file, one can get model parameter estimates and inference results such as eQTL detection based on local false discovery rates, tissue specificity assessment, and 2-dimensional scatter plots of eQTL discoveries. Please find more information in our manuscript at arXiv:1311.2948.

Update 2014.11.3: EM1_XXX is the version with no constraints on the diagonal value of Delta (i.e., commented out the constraint step)