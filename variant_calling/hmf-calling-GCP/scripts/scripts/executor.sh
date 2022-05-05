#!/bin/bash
set -e

cram_path_normal=$1
cram_path_tumoral=$2
bucket=$3
cpus=$4

outfolder=/workspace/analysis/output
mkdir -p ${outfolder}

#gcloud compute instances add-metadata ${CLOUD_INSTANCE} --zone=${CLOUD_ZONE} --metadata STATUSPIPELINE=FIRSTDOWNLOADING
echo "Preprocess" > ${outfolder}/report.txt
uptime >> ${outfolder}/report.txt

bash /workspace/scripts/preprocess.sh ${cram_path_normal} ${cram_path_tumoral}

cram_normal=$(basename ${cram_path_normal})
cram_tumoral=$(basename ${cram_path_tumoral})

bash /workspace/scripts/analysis.sh ${cram_normal} ${cram_tumoral} ${cpus}

echo "Postprocess" >> ${outfolder}/report.txt
uptime >> ${outfolder}/report.txt

outname="${cram_normal}_${cram_tumoral}"
bash /workspace/scripts/postprocess.sh ${bucket} ${outname}
