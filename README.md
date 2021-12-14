# Manipulate and check pedigree file and calculate the triplet format of pedigree relationship matrix (A-matrix)
The R script in this repository evaluates the pedigree data based on the following process:
1) Pedigree Preparation using **prePed** function in **optiSel** package  
2) Checking the completeness of the pedigree for the intended evaluation and identifying individuals with sufficient pedigree information. The following parameters were computed for each individual during the pedigree investigation:
-   **equiGen**: Number of equivalent complete generations, calculated as the sum of the proportion of known ancestors over all generations traced
-   **fullGen**: Number of fully traced generations
-   **maxGen**: Number of maximum generations traced
-   **PCI**: Index of pedigree completeness, which is the harmonic mean of the pedigree completenesses of the parents (MacCluer et al, 1983)
-   **Inbreeding**: Pedigree-based inbreeding coefficient
3) Calculating pedigree-based relationship matrix (A matrix) using **nadiv** package and converting the sparse matrix to triplet format. 
