#!/bin/bash
set -e

type="${1:-scripts}"

# Install bzip2
sudo apt install -y bzip2 git

# Create the appropriate folders
sudo mkdir -p /workspace
sudo chmod go+w /workspace
mkdir -p /workspace/soft
mkdir -p /workspace/analysis
mkdir -p /workspace/data
mkdir -p /workspace/scripts

mv scripts/* /workspace/scripts/

if [ "${type}" == "remote" ]
then
	find /workspace/scripts/ -type f -name '*.sh' | xargs sed -i 's/#gcloud/gcloud/g'
fi

## Download the required data files mentioned in https://github.com/hartwigmedical/gridss-purple-linx
#cd /workspace/data
#wget "https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMFTools-Resources%2FGRIDSS-Purple-Linx-Docker&files=gridss-purple-linx-hg19-refdata-Dec2019.tar.gz" -O gridss-purple-linx-hg19-refdata.tar.gz
#tar -xvf gridss-purple-linx-hg19-refdata.tar.gz hg19/refgenomes/Homo_sapiens.GRCh37.GATK.illumina
## Clean up
#rm gridss-purple-linx-hg19-refdata.tar.gz
# Copy required files
mkdir -p /workspace/data/
gsutil -m cp -r gs://files_bbglab/* /workspace/data/

# Download and install strelka
cd /workspace/soft
wget https://github.com/Illumina/strelka/releases/download/v2.9.2/strelka-2.9.2.centos6_x86_64.tar.bz2
tar xvjf strelka-2.9.2.centos6_x86_64.tar.bz2
# Test strelka
bash strelka-2.9.2.centos6_x86_64/bin/runStrelkaSomaticWorkflowDemo.bash
bash strelka-2.9.2.centos6_x86_64/bin/runStrelkaGermlineWorkflowDemo.bash
# Clean up
rm -r strelkaGermlineDemoAnalysis
rm -r strelkaSomaticDemoAnalysis
rm strelka-2.9.2.centos6_x86_64.tar.bz2

# Get miniconda
cd /workspace/soft
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# Install it and during installation set the installation directory to /workspace/soft/miniconda3 and run conda init
bash Miniconda3-latest-Linux-x86_64.sh -p /workspace/soft/miniconda3 -b

# Create a source file for conda init
cat > /workspace/scripts/condainit << EOF

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="\$('/workspace/soft/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ \$? -eq 0 ]; then
    eval "\$__conda_setup"
else
    if [ -f "/workspace/soft/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/workspace/soft/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/workspace/soft/miniconda3/bin:\$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

EOF

rm Miniconda3-latest-Linux-x86_64.sh

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed
chmod +x bigWigAverageOverBed
git clone https://github.com/parklab/MosaicForecast.git

# Install conda packages
source /workspace/scripts/condainit
conda create -n samtools -y -c bioconda -c conda-forge samtools

conda create -n tabix -y -c bioconda -c conda-forge tabix

conda create -n bedtools -y -c bioconda -c conda-forge bedtools

conda env create --name MF --file MosaicForecast/environment.yaml