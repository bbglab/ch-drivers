# HMF and TCGA processing

After the variant calling is performed, the data must be postprocessed. We have different scripts for the different sequencing projects we are analysing.
There are harcoded paths in each script that need to be changed to your current setup. 
Packages needed to run these scripts include pandas, numpy, tqdm, pybedtools, bgreference, pyliftover, bgparsers


Files needed:

- gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz (Download from ENSEMBl)
- gnomad.genomes.r2.1.1.sites.vcf.bgz  (Download from ENSEMBl)
- MuTect2.PON.4136.vcf.gz (Download from https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)
- Homo_sapiens_assembly38.sdust.30.gz or Homo_sapiens.GRCh37.GATK.illumina.fasta.sdust.30.gz (Run https://pubmed.ncbi.nlm.nih.gov/16796549/ on fasta files)
- k36.umap.bed.noheader.gz  (Download from https://bismap.hoffmanlab.org/, remove the header)
- SageGermlinePon.hg19.vcf.gz (Download from HMF Tools)

HMF
`
python hmf_process.py $FILE $OUTPATH 
`
TCGA 
`
python tcga_process.py $FILE $OUTPATH
`

Then merge the data (changing the corresponding paths) with 

`python merging.py `


We will then apply filters on these merged files. These filters include 

`
python process_mutations.py HMF_rev_muts_gnomad.hg19.tsv.gz HMF vaf_filter/ 0
`
`
python process_mutations.py TCGA_rev_muts_gnomad.hg38.tsv.gz TCGA vaf_filter/ 0
`

In order to run the script, several files are needed, including

- snp151Common.hg19.vep.txt.sorted.bgz from dbSNP
- 00-common_all.grch37.vcf.gz from dbSNP
- metadata.tsv (from HMF when you request access)
- clinical_PANCAN_patient_with_followup_151210.tsv (from https://gdc.cancer.gov/about-data/publications/pancanatlas ) 
- repeat_masker_hg38.gz (from UCSC)
- segmental_duplications_hg38.gz (from UCSC)
- simple_repeats_hg38.gz (from UCSC)

This will generate the input files for INTOGEN and other files needed for posterior analyses.


We then run vep. First we generate the correct VEP input

` 
python vep_annotate.py
`

Then we run VEP (we are using Singularity, but there are other options including conda)

`
./run_vep.sh
`

Finally, we will postprocess VEP output:

`
python vep_postprocess.py
` 

# MSKCC impact

We will process the data obtained from (https://www.cbioportal.org/study/summary?id=msk_ch_2020)

`
python impact_preformat.py
`