 as always change your paths to your own singularity image (you can download it from ensemble) and your own vep cache

# VEP container can be downloaded from https://hub.docker.com/r/ensemblorg/ensembl-vep
singularity exec /workspace/datasets/intogen/oriolRun/containers/vep.simg vep -i /workspace/projects/reverse_calling/rebuttle/final_run/HMF_full.vep.sorted.tsv.gz \
-o STDOUT --assembly GRCh37 --no_stats --cache --offline --symbol --protein --tab --canonical --dir /workspace/datasets/vep/ --fork 32\
|  grep -v ^## |  gzip > /workspace/projects/reverse_calling/rebuttle/final_run/HMF_full.vep.sorted.tsv.out.gz


singularity exec /workspace/datasets/intogen/oriolRun/containers/vep.simg vep -i /workspace/projects/reverse_calling/rebuttle/final_run/TCGA_full.vep.sorted.tsv.gz \
-o STDOUT --assembly GRCh38 --no_stats --cache --offline --symbol --protein --tab --canonical --dir /workspace/datasets/intogen/oriolRun/datasets/vep --fork 32\
|  grep -v ^## | gzip > /workspace/projects/reverse_calling/rebuttle/final_run/TCGA_full.vep.sorted.tsv.out.gz