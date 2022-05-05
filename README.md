# Discovering the drivers of clonal hematopoiesis

This is the code to reproduce the analyses described in **Discovering the drivers of clonal hematopoiesis**. If you use the code and the data please cite:

   Oriol Pich, Iker Reyes-Salazar, Abel Gonzalez-Perez, Nuria Lopez-Bigas, **Discovering the drivers of clonal hematopoiesis**, Nature Communications (2022)


## How to access the data
It is important to remark that in order to reproduce the analyses several access-controlled files are needed.

The sequencing data to carry out the reverse calling of blood somatic mutations (and germline variants across donors) is available via dbGaP (TCGA; phs000178.v11.p8; https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000178.v11.p8) and HMF (https://hartwigmedical.github.io/documentation/data-access-request-application.html, version DR110).  Access to these protected data must be requested from TCGA and HMF. The procedure and conditions to access these datasets are detailed in the sites referenced above. Gene expression in bone marrow CD34+ cells are available at The Gene Expression Omnibus (GSE96811; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96811). H3K27ac ChIP data for CD34+ samples are available from ENCODE (https://www.encodeproject.org/). Mutations in CH drivers across hematopoietic malignancies are available from IntOGen (intogen.org). Disease-related variants are available from ClinVar (https://ftp.ncbi.nlm.nih.gov/pub/clinvar/).

We have prepared flat files containing the set of blood somatic mutations identified in both datasets and have made them available through HMF and dbGaP following the same procedure to access the original datasets. HMF blood somatic mutations are available as part of the data access request to HMF (see above). TCGA blood somatic mutations are available through dbGaP (phs002867.v1.p1; https://www.ncbi.nlm.nih.gov/projects/gapprev/gap/cgi-bin/study.cgi?study_id=phs002867.v1.p1) to researchers who have obtained permission to access protected TCGA data. Panel-sequenced data from the IMPACT targeted cohort is available through cBioPortal (https://www.cbioportal.org/study/summary?id=msk_ch_2020). Other datasets employed in specific analyses are described in prior sections of these Methods and within the code repository. 


## How to reproduce the analyses
In order to reproduce the analyses from scratch, first start running the scripts located in the **variant_calling** folder, and then continue in the **src** folder.

The code to reproduce all analyses and figures in the paper is in the **figures** folder.
