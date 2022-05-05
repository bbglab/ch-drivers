#!/bin/bash
set -e


cram_normal=$1
cram_tumoral=$2
folder="${3:-/workspace/analysis/input}"

mkdir -p ${folder}

gsutil -u ${CLOUD_PROJECT} \
	cp ${cram_normal} ${folder}

#gcloud compute instances add-metadata ${CLOUD_INSTANCE} --zone=${CLOUD_ZONE} --metadata STATUSPIPELINE=DOWNLOADED1

gsutil -u ${CLOUD_PROJECT} \
	cp ${cram_tumoral} ${folder}

#gcloud compute instances add-metadata ${CLOUD_INSTANCE} --zone=${CLOUD_ZONE} --metadata STATUSPIPELINE=DOWNLOADED

gsutil -u ${CLOUD_PROJECT} cp ${cram_normal}.crai ${folder}
gsutil -u ${CLOUD_PROJECT} cp ${cram_tumoral}.crai ${folder}

#gcloud compute instances add-metadata ${CLOUD_INSTANCE} --zone=${CLOUD_ZONE} --metadata STATUSPIPELINE=DOWNLOADED