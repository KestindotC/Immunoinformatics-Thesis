Immunoinformatics Thesis Concept Map
=============


Characterizing potential tumor immune association in colorectal carcinoma(CRC).
Base on immune cell infiltration among the colorectal cancer cohorts, 
we aim to find out which immune characteristics were associated with the abundance of TILs.


Immunoinformatics analysit were dedicated in exploring some features of tumor immnity.   
We took consider in following issues, and exploring association between features based on the TILs expression:
1. HLA alleles
2. Neoantigen
3. Gene mutation (SNV only)


## Tumor-infiltrating lymphocytes
Following RNA-sequencing data analysis pipeline, we calculated gene set score of many immune subtypes.
**Figure 1** shows clustering of total 623 CRC samples.



## HLA typing
**POLYSOLVER** were applied to obtain HLA alleles in each patient.



## Neoantigen Prediction
netMHCpan3.0 [link] were applied in our scripts, Neopack.py.
filtered.txt files were filtered by strong binding or weak binding.

