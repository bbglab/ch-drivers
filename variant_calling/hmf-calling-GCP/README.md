# Prepare VM on GCP to run the Hartwig analysis

## Prepare the base image

Create a VM with very little resources based on Debian.

Remember to set up the zone and the appropriate service account
(although it might not be used in this phase).

Copy the scripts folder to the running instance:

```bash
gcloud compute scp --recurse scripts <user>@<instance-name>:~
```

Execute the install script in the VM to check for errors.
Two install options are available: local and remote.
Local is meant for use of the as a user interactively:

```bash
cd scripts
bash install.sh
```

Remote is for make use of the service account and
the metadata to retrieve the status of the VM:

```bash
cd scripts
bash install.sh remote
```

Now you can create an image from this VM.

## Executing a run as user

To login in the VM with your service account use

```bash
gcloud auth login
```

and copy the verification code

You can execute a simple run using the executor.sh script as:

```bash
bash /workspace/scripts/executor.sh <normal-bam-path> <tumoral-bam-path> <your-bucket> <cpus>
```

The ``CLOUD_PROJECT`` environment variable must be defined to be able
to copy the bam files.

## Executing using gcloud API

To use the API you need: click and google-api-python-client

```bash
pip install --upgrade google-api-python-client
```

## Gnomad

Download gnomad from: ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/
Files needed are: ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz
and ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz.tbi

```bash
bcftools norm -m-snps gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz |\ 
    bcftools view  --min-af 0.01 --types snps | grep -v "^#" |\
    awk -F $'\t' 'BEGIN {OFS = FS} {print $1,$2,substr($4,1,1),substr($5,1,1)}' |\ 
    sort -k1V -k2n |\ 
    bgzip > gnomad.genomes.r2.0.1.site.noVEP.0.01.tsv.bgz
tabix -s 1 -b 2 -e 2 gnomad.genomes.r2.0.1.site.noVEP.0.01.tsv.bgz
```