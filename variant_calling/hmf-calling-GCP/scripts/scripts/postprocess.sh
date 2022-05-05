#!/bin/bash
set -e

bucket=$1
outname=$2
folder=${3:-/workspace/analysis/output}

# compress files
gzip ${folder}/filtered_vars.tsv
gzip ${folder}/mosaic/*

# copy results to bucket
gsutil -m cp -r ${folder} $bucket/$outname/

# clean all

#gcloud compute instances add-metadata ${CLOUD_INSTANCE} --zone=${CLOUD_ZONE} --metadata STATUSPIPELINE=DONE