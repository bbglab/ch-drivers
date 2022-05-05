#!/bin/bash
set -e

source /workspace/scripts/condainit
conda activate samtools

cram_normal=$1
cram_tumoral=$2
cpus=$3

infolder=${4:-/workspace/analysis/input}
outfolder=${5:-/workspace/analysis/output}

tmpdir=/workspace/analysis/tmp
mkdir -p ${tmpdir}

REF=/workspace/data/hg19/refgenomes/Homo_sapiens.GRCh37.GATK.illumina.fasta
GNOMAD=/workspace/data/gnomad.genomes.r2.0.1.site.noVEP.0.01.tsv.bgz
K24=/workspace/data/k24.umap.wg.bw

touch ${GNOMAD}.tbi

mkdir -p ${outfolder}

file_normal=${cram_normal}
file_tumoral=${cram_tumoral}

#----------------
# Variant calling
#----------------

echo "Strelka" >> ${outfolder}/report.txt
uptime >> ${outfolder}/report.txt

sample_normal=`samtools view -H ${infolder}/${file_normal} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq`
sample_tumoral=`samtools view -H ${infolder}/${file_tumoral} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq`


# TODO remove the last subfolder (check below also)
out_strelka="${tmpdir}/strelka/${sample_tumoral}_${sample_normal}"
mkdir -p ${out_strelka}

python2 /workspace/soft/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
	--normalBam ${infolder}/${file_normal} --tumorBam ${infolder}/${file_tumoral} \
	--referenceFasta ${REF} \
	--runDir ${out_strelka}

python2 ${out_strelka}/runWorkflow.py -m local -j ${cpus}

strelka_results="${outfolder}/strelka"
mkdir -p ${strelka_results}

cp -r ${out_strelka}/results ${strelka_results}/
cp ${out_strelka}/*.txt ${strelka_results}/

echo "STRELKA done"


#----------------------
# Filter strelka output
#----------------------
conda deactivate
conda activate tabix

echo "Strelka vars" >> ${outfolder}/report.txt
uptime >> ${outfolder}/report.txt

strelka_vars=${tmpdir}/strelka_vars.tsv

zcat ${out_strelka}/results/variants/somatic.snvs.vcf.gz |\
 	grep -v "#" | grep -w PASS | cut -f1,2,4,5 \
 	> ${strelka_vars}

echo `wc -l ${strelka_vars}` " variants"

splits_folder=${tmpdir}/splits
mkdir -p ${splits_folder}

divisor=$((${cpus}-1))

split -l$((`wc -l < ${strelka_vars}`/${divisor})) ${strelka_vars} \
	${splits_folder}/strelka_vars.split. -d

splits_filtered_folder=${tmpdir}/splits_filtered
mkdir -p ${splits_filtered_folder}

for f in ${splits_folder}/*.split*
do
	name=${f##*.}
	awk -v tabixfile=${GNOMAD} -f /workspace/scripts/filter.awk $f > ${splits_filtered_folder}/${name}.tsv &
done

wait

cat ${splits_filtered_folder}/*.tsv > ${outfolder}/filtered_vars.tsv
echo "# variants" >> ${outfolder}/report.txt
echo `wc -l ${outfolder}/filtered_vars.tsv` >> ${outfolder}/report.txt

awk -F $'\t' 'BEGIN {OFS = FS} {print $1,$2-1,$2}' ${outfolder}/filtered_vars.tsv \
	> ${tmpdir}/filtered_vars.bed

conda deactivate
conda activate bedtools

bedtools merge -d 152 -i ${tmpdir}/filtered_vars.bed |\
	awk '{print $1":"$2"-"$3}' > ${tmpdir}/vars4download.txt

splits_download=${tmpdir}/vars4download_splits
mkdir -p ${splits_download}
split -l$((`wc -l < ${tmpdir}/vars4download.txt`/${divisor})) ${tmpdir}/vars4download.txt \
	${splits_download}/split. -d

awk -F $'\t' 'BEGIN {OFS = FS} {print $1,$2-2002,$2+2001}' ${outfolder}/filtered_vars.tsv \
	> ${tmpdir}/filtered_vars_expanded.bed

bedtools merge -d 152 -i ${tmpdir}/filtered_vars_expanded.bed |\
	awk '{print $1":"$2"-"$3}' > ${tmpdir}/vars4phase.txt

splits_phase=${tmpdir}/vars4phase_splits
mkdir -p ${splits_phase}
split -l$((`wc -l < ${tmpdir}/vars4phase.txt`/${divisor})) ${tmpdir}/vars4phase.txt \
	${splits_phase}/split. -d

#---------------------------
# Reads supporting variants
#---------------------------
conda deactivate
conda activate samtools

echo "Support reads" >> ${outfolder}/report.txt
uptime >> ${outfolder}/report.txt

variants_healthy_folder=${tmpdir}/healthy
mkdir -p ${variants_healthy_folder}

for f in ${splits_download}/*
do
	name=${f##*.}
	cat $f | xargs samtools view -b -T ${REF} -o ${variants_healthy_folder}/${name}.bam ${infolder}/${file_normal}  &
done

wait

samtools merge -@ ${cpus} ${outfolder}/variants_healthy_reads.bam ${variants_healthy_folder}/*

variants_tumor_folder=${tmpdir}/tumoral
mkdir -p ${variants_tumor_folder}

for f in ${splits_download}/*
do
	name=${f##*.}
	cat $f | xargs samtools view -b -T ${REF} -o ${variants_tumor_folder}/${name}.bam ${infolder}/${file_tumoral} &
done

wait

samtools merge -@ ${cpus} ${outfolder}/variants_tumor_reads.bam ${variants_tumor_folder}/*

variants_phase_folder=${tmpdir}/phase
mkdir -p ${variants_phase_folder}

for f in ${splits_phase}/*
do
	name=${f##*.}
	cat $f | xargs xargs samtools view -b -T ${REF} -o ${variants_phase_folder}/${name}.bam ${infolder}/${file_tumoral} &
done

wait

samtools merge -@ ${cpus} ${tmpdir}/phase_unsorted.bam ${variants_phase_folder}/*
samtools sort -@ ${cpus} -o ${tmpdir}/phase.bam ${tmpdir}/phase_unsorted.bam
samtools index -@ ${cpus} -b ${tmpdir}/phase.bam

#-----------------------------------
# add reads mapping CH related genes
#-----------------------------------

#-----------------
# Run Mosaic
#-----------------
#gcloud compute instances add-metadata ${CLOUD_INSTANCE} --zone=${CLOUD_ZONE} --metadata STATUSPIPELINE=RUNNING_MOSAIC
conda deactivate
# allow bigWigAverageOverBed to be called directly
export PATH="/workspace/soft:$PATH"
conda activate MF
out_mosaic="${tmpdir}/mosaic"
mkdir -p ${out_mosaic}

echo "Mosaic phase" >> ${outfolder}/report.txt
uptime >> ${outfolder}/report.txt

mosaic_vars_input=${tmpdir}/mosaic_input
awk -v sample=phase -F $'\t' 'BEGIN {OFS = FS} {print $1,$2-1,$2,$3,$4,sample}' \
	${outfolder}/filtered_vars.tsv > ${mosaic_vars_input}

python /workspace/soft/MosaicForecast/Phase.py \
	${tmpdir} ${out_mosaic} ${REF} \
	${mosaic_vars_input} 20 \
	${K24} ${cpus} bam

echo "Mosaic features" >> ${outfolder}/report.txt
uptime >> ${outfolder}/report.txt

python /workspace/soft/MosaicForecast/ReadLevel_Features_extraction.py \
	${mosaic_vars_input} ${out_mosaic}/features.txt  \
	${tmpdir}  ${REF} \
	${K24} ${cpus} bam


echo "Mosaic predict" >> ${outfolder}/report.txt
uptime >> ${outfolder}/report.txt

Rscript /workspace/soft/MosaicForecast/Prediction.R ${out_mosaic}/features.txt \
	/workspace/soft/MosaicForecast/models_trained/50xRFmodel_addRMSK_Refine.rds \
	Refine ${out_mosaic}/predictions_test.tsv

mosaic_results="${outfolder}/mosaic"
mkdir -p ${mosaic_results}

for f in `find ${out_mosaic} -maxdepth 1 -type f -not -name "*.tmp"`
do
	name=$(basename $f)
	cp $f ${mosaic_results}/${name}
done
