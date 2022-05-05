TUMORAL=$(curl http://metadata/computeMetadata/v1/instance/attributes/tumoral -H "Metadata-Flavor: Google")
NORMAL=$(curl http://metadata/computeMetadata/v1/instance/attributes/normal -H "Metadata-Flavor: Google")

CS_BUCKET=$(curl http://metadata/computeMetadata/v1/instance/attributes/bucket -H "Metadata-Flavor: Google")
CPUS=$(curl http://metadata/computeMetadata/v1/instance/attributes/cpus -H "Metadata-Flavor: Google")

project=$(curl -X GET http://metadata/computeMetadata/v1/instance/attributes/project -H 'Metadata-Flavor: Google')
instance=$(curl -X GET http://metadata/computeMetadata/v1/instance/attributes/name -H 'Metadata-Flavor: Google')
zone=$(curl -X GET http://metadata/computeMetadata/v1/instance/attributes/zone -H 'Metadata-Flavor: Google')

export CLOUD_PROJECT=${project}
export CLOUD_INSTANCE=${instance}
export CLOUD_ZONE=${zone}

bash /workspace/scripts/executor.sh $NORMAL $TUMORAL $CS_BUCKET $CPUS

if [ $? -eq 0 ]; then
    sleep 5s
    gcloud compute instances add-metadata ${CLOUD_INSTANCE} --zone=${CLOUD_ZONE} --metadata finished_task=succeed
else
    sleep 5s
    gcloud compute instances add-metadata ${CLOUD_INSTANCE} --zone=${CLOUD_ZONE} --metadata finished_task=FAIL
fi

# sleep so that we can check the metadata