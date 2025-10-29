#!/usr/bin/env bash

RESOURCES_DIR=${1:-"Resources"};
BEDTOOLS_DOCKER=${2:-"levanonlab/bedtools:2.31.1"}
BIGBED_TO_BED_DOCKER=${3:-"levanonlab/bigbedtobed:482--h0b57e2e_0"}
RNA_EDITING_INDEX_DOCKER=${4:-'levanonlab/rna-editing-index-lite:1.0'}
DBS_PATHS_INI=${5:-"${RESOURCES_DIR}/ResourcesPaths.ini"};
DONT_DOWNLOAD=${6:-false}
DONT_WRITE=${7:-false}
DONT_GENERATE_GENOME_INDEXES=${8:-false}

#---------------------------------------------------------------------------
# Constants
#---------------------------------------------------------------------------
HUMAN="HomoSapiens"
MURINE="MusMusculus"

HG38_FTP_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/"
HG38_FTP_GENOME_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/"
HG38_FTP_GENES_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/"

MM10_FTP_URL="http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/"
MM10_FTP_GENOME_URL="http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/"
MM10_FTP_GENES_URL="http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/"

MM10_GENE_EXPRESSION_FTP="https://hgdownload.soe.ucsc.edu/gbdb/mm10/tabulamuris/"
MM10_GENES_EXPRESSION_TABLE_FILE="barChart.bb"

GENOME_DIR="Genomes"
HUMAN_GENOME_DIR="${RESOURCES_DIR}/${GENOME_DIR}/${HUMAN}"
HG38_GENOME_FASTA_FILE="hg38.fa.gz"
HG38_GENOME_FASTA="ucscHg38Genome.fa"

MURINE_GENOME_DIR="${RESOURCES_DIR}/${GENOME_DIR}/${MURINE}"
MM10_GENOME_FASTA_FILE="mm10.fa.gz" 
MM10_GENOME_FASTA="ucscMm10Genome.fa"


TRANSCRIPTOME_DIR="Transcriptomes"
HUMAN_TRANSCRIPTOME_DIR="${RESOURCES_DIR}/${TRANSCRIPTOME_DIR}/${HUMAN}"
HG38_TRANSCRIPTOME_FASTA_FILE="refMrna.fa.gz"
HG38_TRANSCRIPTOME_FASTA="ucscHg38RefSeqTranscriptome.fa"

MURINE_TRANSCRIPTOME_DIR="${RESOURCES_DIR}/${TRANSCRIPTOME_DIR}/${MURINE}"
MM10_TRANSCRIPTOME_FASTA_FILE="refMrna.fa.gz" 
MM10_TRANSCRIPTOME_FASTA="ucscMm10RefSeqTranscriptome.fa"


REGIONS_DIR="Regions"
HUMAN_REGIONS_DIR="${RESOURCES_DIR}/${REGIONS_DIR}/${HUMAN}"
HG38_REGIONS_FILE="ucscHg38Alu.bed.gz"
HG38_REGIONS_TABLE_FILE="rmsk.txt.gz"

MURINE_REGIONS_DIR="${RESOURCES_DIR}/${REGIONS_DIR}/${MURINE}"
MM10_REGIONS_FILE="ucscMM10SINE_B1_B2.bed.gz"
MM10_REGIONS_TABLE_FILE="rmsk.txt.gz"


SNPS_DIR="SNPs"
HUMAN_SNPS_DIR="${RESOURCES_DIR}/${SNPS_DIR}/${HUMAN}"
HG38_SNPS_FILE="ucscHg38CommonGenomicSNPs151.bed.gz"
HG38_SNPS_TABLE_FILE="snp151Common.txt.gz"

MURINE_SNPS_DIR="${RESOURCES_DIR}/${SNPS_DIR}/${MURINE}"
MM10_SNPS_FILE="ucscMM10CommonGenomicSNPs142.bed.gz"
MM10_SNPS_TABLE_FILE="snp142Common.txt.gz"


REFSEQ_DIR="RefSeqAnnotations"
HUMAN_REFSEQ_DIR="${RESOURCES_DIR}/${REFSEQ_DIR}/${HUMAN}"
HG38_REFSEQ_TABLE_FILE="ncbiRefSeqCurated.txt.gz"
HG38_REFSEQ_FILE="ucscHg38RefSeqCurated.bed.gz"
HG38_REFSEQ_FILE_FINAL="ucscHg38RefSeqCurated.tsv"
HG38_REFSEQ_GTF_ZIP_FILE="hg38.ncbiRefSeq.gtf.gz"
HG38_REFSEQ_GTF="hg38.ncbiRefSeq.gtf"
HG38_REFSEQ_GTF_FINAL="hg38.ncbiRefSeq.noTranscriptVersion.gtf"

MURINE_REFSEQ_DIR="${RESOURCES_DIR}/${REFSEQ_DIR}/${MURINE}"
MM10_REFSEQ_TABLE_FILE="ncbiRefSeqCurated.txt.gz"
MM10_REFSEQ_FILE="ucscMM10RefSeqCurated.bed.gz"
MM10_REFSEQ_FILE_FINAL="ucscMM10RefSeqCurated.tsv"
MM10_REFSEQ_GTF_ZIP_FILE="mm10.ncbiRefSeq.gtf.gz"
MM10_REFSEQ_GTF_FILE="mm10.ncbiRefSeq.gtf"

GENES_EXPRESSION_DIR="GenesExpression"
HUMAN_GENES_EXPRESSION_DIR="${RESOURCES_DIR}/${GENES_EXPRESSION_DIR}/${HUMAN}"
HG38_GENES_EXPRESSION_FILE="ucscHg38GTExGeneExpression.bed.gz"
HG38_GENES_EXPRESSION_TABLE_FILE="gtexGene.txt.gz"

MURINE_GENES_EXPRESSION_DIR="${RESOURCES_DIR}/${GENES_EXPRESSION_DIR}/${MURINE}"
MM10_GENES_EXPRESSION_FILE="ucscMM10GeneExpression.bed.gz"

#---------------------------------------------------------------------------
# Clone needed scripts
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# Downloading
#---------------------------------------------------------------------------

if [ "${DONT_DOWNLOAD}" = false ]
then

  # clean folders from previous runs
  echo "Info: Trying to delete old resource files"
  find "${RESOURCES_DIR}" -type f -delete

  mkdir -p "${HUMAN_GENOME_DIR}"
  mkdir -p "${MURINE_GENOME_DIR}"

  mkdir -p "${HUMAN_REGIONS_DIR}"
  mkdir -p "${MURINE_REGIONS_DIR}"

  mkdir -p "${HUMAN_SNPS_DIR}"
  mkdir -p "${MURINE_SNPS_DIR}"

  mkdir -p "${HUMAN_REFSEQ_DIR}"
  mkdir -p "${MURINE_REFSEQ_DIR}"

  mkdir -p "${HUMAN_GENES_EXPRESSION_DIR}"
  mkdir -p "${MURINE_GENES_EXPRESSION_DIR}"


  echo "Info: Started Downloading UCSC Resources."

  #---------------------------------------------------------------------------
  # HG38
  #---------------------------------------------------------------------------
  echo "Info: Started Downloading Hg38 Files:"

  # Genome
  echo "Downloading Hg38 Genome: ${HG38_FTP_GENOME_URL}${HG38_GENOME_FASTA_FILE}"
  wget "${HG38_FTP_GENOME_URL}${HG38_GENOME_FASTA_FILE}"  --directory-prefix="${HUMAN_GENOME_DIR}"
  echo "Saving Hg38 Genome Under: ${HUMAN_GENOME_DIR}/${HG38_GENOME_FASTA}"
  gunzip -c "${HUMAN_GENOME_DIR}/${HG38_GENOME_FASTA_FILE}" > "${HUMAN_GENOME_DIR}/${HG38_GENOME_FASTA}"
  rm "${HUMAN_GENOME_DIR}/${HG38_GENOME_FASTA_FILE}"
  echo "Done Processing Hg38 Genome"

  # Transcriptome
  echo "Downloading Hg38 Transcriptome: ${HG38_FTP_GENOME_URL}${HG38_TRANSCRIPTOME_FASTA_FILE}"
  wget "${HG38_FTP_GENOME_URL}${HG38_TRANSCRIPTOME_FASTA_FILE}"  --directory-prefix="${HUMAN_TRANSCRIPTOME_DIR}"
  echo "Saving Hg38 Transcriptome Under: ${HUMAN_TRANSCRIPTOME_DIR}/${HG38_TRANSCRIPTOME_FASTA}"
  gunzip -c "${HUMAN_TRANSCRIPTOME_DIR}/${HG38_TRANSCRIPTOME_FASTA_FILE}" > "${HUMAN_TRANSCRIPTOME_DIR}/${HG38_TRANSCRIPTOME_FASTA}"
  rm "${HUMAN_TRANSCRIPTOME_DIR}/${HG38_TRANSCRIPTOME_FASTA_FILE}"
  echo "Done Processing Hg38 Transcriptome"
  
  # Repeats Regions
  echo "Downloading Hg38 Alu Repeats Table ${HG38_FTP_URL}${HG38_REGIONS_TABLE_FILE}"
  wget "${HG38_FTP_URL}${HG38_REGIONS_TABLE_FILE}"  --directory-prefix="${HUMAN_REGIONS_DIR}"
  echo "Processing Hg38  Alu Repeats Table ${HG38_REGIONS_TABLE_FILE}"
  zcat "${HUMAN_REGIONS_DIR}/${HG38_REGIONS_TABLE_FILE}"| awk '{OFS ="\t"}($13 ~/Alu/ && $6 !~/_/) {print $6,$7,$8}' > "${HUMAN_REGIONS_DIR}/${HG38_REGIONS_TABLE_FILE}.tmp"
  docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v "$(pwd)/${HUMAN_REGIONS_DIR}":/data -t ${BEDTOOLS_DOCKER} /bin/bash -c "bedtools sort -i /data/${HG38_REGIONS_TABLE_FILE}.tmp | bedtools merge -i stdin | gzip > /data/${HG38_REGIONS_FILE}"
  rm "${HUMAN_REGIONS_DIR}/${HG38_REGIONS_TABLE_FILE}" "${HUMAN_REGIONS_DIR}/${HG38_REGIONS_TABLE_FILE}.tmp"
  echo "Done Processing Hg38 Alu Repeats Table ${HG38_REGIONS_TABLE_FILE}"
    # Generate Indexes
  if [ "${DONT_GENERATE_GENOME_INDEXES}" = false ]
  then
    # JSD path is only genome name to match our Nextflow workflow requisits
    echo "Attempting to Create Genome Index of Hg38 Alu Repeats ${HG38_REGIONS_FILE}"
    CMD="cd / ; /usr/bin/java -jar /bin/AEI/RNAEditingIndexer/lib/EditingIndexJavaUtils.jar GenerateIndex -i ${HUMAN_REGIONS_DIR}/${HG38_REGIONS_FILE} -g ${HUMAN_GENOME_DIR}/${HG38_GENOME_FASTA} -o ${HUMAN_REGIONS_DIR}/${HG38_GENOME_FASTA}.${HG38_REGIONS_FILE}.GenomeIndex.jsd -b bedtools"
    # CMD="pwd"
    docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v "$(pwd)/${RESOURCES_DIR}":"/${RESOURCES_DIR}" -t ${RNA_EDITING_INDEX_DOCKER} /bin/bash -c "${CMD}"
    echo "Done Creating Genome Index of Hg38 Alu Repeats ${HUMAN_REGIONS_DIR}/${HG38_GENOME_FASTA}.${HG38_REGIONS_FILE}.GenomeIndex.jsd"
  fi
  
  # SNPs
  echo "Downloading Hg38 Common Genomic SNPs Table ${HG38_FTP_URL}${HG38_SNPS_TABLE_FILE}"
  wget "${HG38_FTP_URL}${HG38_SNPS_TABLE_FILE}"  --directory-prefix="${HUMAN_SNPS_DIR}"
  echo "Processing Hg38Common Genomic SNPs Table ${HG38_SNPS_TABLE_FILE}"
  zcat "${HUMAN_SNPS_DIR}/${HG38_SNPS_TABLE_FILE}" | awk '{OFS ="\t"}($11=="genomic") {print $2,$3,$4,$7,$9,$10,$16,$25}'| gzip > "${HUMAN_SNPS_DIR}/${HG38_SNPS_FILE}"
  rm "${HUMAN_SNPS_DIR}/${HG38_SNPS_TABLE_FILE}"
  echo "Done Processing Hg38 Common Genomic SNPs Table ${HG38_SNPS_TABLE_FILE}"

  # RefSeq Annotations
  echo "Downloading Hg38 RefSeq Curated Table ${HG38_FTP_URL}${HG38_REFSEQ_TABLE_FILE}"
  wget "${HG38_FTP_URL}${HG38_REFSEQ_TABLE_FILE}"  --directory-prefix="${HUMAN_REFSEQ_DIR}"
  echo "Processing Hg38 RefSeq Curated Table ${HG38_REFSEQ_TABLE_FILE}"
  zcat "${HUMAN_REFSEQ_DIR}/${HG38_REFSEQ_TABLE_FILE}"| awk '{OFS ="\t"} {print $3,$5,$6,$2,$13,$4,$10,$11}' |gzip > "${HUMAN_REFSEQ_DIR}/${HG38_REFSEQ_FILE}"
  echo "Generating Hg38 RefSeq Curated Transcript-to-Gene file from ${HG38_FTP_URL}${HG38_REFSEQ_FILE}"
  zcat "${HUMAN_REFSEQ_DIR}/${HG38_REFSEQ_FILE}" | awk 'BEGIN {FS=OFS="\t"; print "Transcript\tGene"} { gsub(/\.[0-9_]+/,"",$4) ; print $4,$5}' > "${HUMAN_REFSEQ_DIR}/${HG38_REFSEQ_FILE_FINAL}"
  rm "${HUMAN_REFSEQ_DIR}/${HG38_REFSEQ_TABLE_FILE}"
  echo "Done Processing Hg38 RefSeq Curated Table ${HG38_REFSEQ_TABLE_FILE}"

  echo "Downloading Hg38 RefSeq GTF Annotation File ${HG38_FTP_GENES_URL}${HG38_REFSEQ_GTF_ZIP_FILE}"
  curl "${HG38_FTP_GENES_URL}${HG38_REFSEQ_GTF_ZIP_FILE}" -o "${HUMAN_REFSEQ_DIR}/${HG38_REFSEQ_GTF_ZIP_FILE}"
  echo "Saving Hg38 RefSeq Table Under: ${HUMAN_REFSEQ_DIR}/${HG38_REFSEQ_GTF_FILE}"
  gunzip "${HUMAN_REFSEQ_DIR}/${HG38_REFSEQ_GTF_ZIP_FILE}"
  awk 'BEGIN {FS=OFS="\t"}{ if ($3 == "transcript") { gsub(/\.[0-9_]+/,"",$9); print } }' "${HUMAN_REFSEQ_DIR}/${HG38_REFSEQ_GTF}" > "${HUMAN_REFSEQ_DIR}/${HG38_REFSEQ_GTF_FINAL}"
  echo "Done Processing Hg38 RefSeq Table ${HG38_REFSEQ_GTF_ZIP_FILE}"

  # Genes Expression
  echo "Downloading Hg38 Genes Expression Table ${HG38_FTP_URL}${HG38_GENES_EXPRESSION_TABLE_FILE}"
  wget "${HG38_FTP_URL}${HG38_GENES_EXPRESSION_TABLE_FILE}"  --directory-prefix="${HUMAN_GENES_EXPRESSION_DIR}"
  echo "Processing Hg38 Genes Expression Table ${HG38_GENES_EXPRESSION_TABLE_FILE}"
  zcat "${HUMAN_GENES_EXPRESSION_DIR}/${HG38_GENES_EXPRESSION_TABLE_FILE}" | awk '{OFS ="\t"} {print $1,$2,$3,$4,$10,$6}'| gzip > "${HUMAN_GENES_EXPRESSION_DIR}/${HG38_GENES_EXPRESSION_FILE}"
  rm "${HUMAN_GENES_EXPRESSION_DIR}/${HG38_GENES_EXPRESSION_TABLE_FILE}"
  echo "Done Processing Hg38 Genes Expression Table ${HG38_GENES_EXPRESSION_TABLE_FILE}"


  #---------------------------------------------------------------------------
  # MM10
  #---------------------------------------------------------------------------
  echo "Info: Started Downloading MM10 Files:"

  # Genome
  echo "Downloading MM10 Genome: ${MM10_FTP_GENOME_URL}${MM10_GENOME_FASTA_FILE}"
  wget "${MM10_FTP_GENOME_URL}${MM10_GENOME_FASTA_FILE}"  --directory-prefix="${MURINE_GENOME_DIR}"
  echo "Saving MM10 Genome Under: ${MURINE_GENOME_DIR}/${MM10_GENOME_FASTA}"
  #tar -xOzf "${MURINE_GENOME_DIR}/${MM10_GENOME_FASTA_FILE}" | cat > "${MURINE_GENOME_DIR}/${MM10_GENOME_FASTA}" # old using TAR
  gunzip -c "${MURINE_GENOME_DIR}/${MM10_GENOME_FASTA_FILE}" > "${MURINE_GENOME_DIR}/${MM10_GENOME_FASTA}"
  rm "${MURINE_GENOME_DIR}/${MM10_GENOME_FASTA_FILE}"
  echo "Done Processing MM10 Genome"

  # Transcriptome
  echo "Downloading Mm10 Transcriptome: ${MM10_FTP_GENOME_URL}${MM10_TRANSCRIPTOME_FASTA_FILE}"
  wget "${MM10_FTP_GENOME_URL}${MM10_TRANSCRIPTOME_FASTA_FILE}"  --directory-prefix="${MURINE_TRANSCRIPTOME_DIR}"
  echo "Saving Mm10 Transcriptome Under: ${MURINE_TRANSCRIPTOME_DIR}/${MM10_TRANSCRIPTOME_FASTA}"
  gunzip -c "${MURINE_TRANSCRIPTOME_DIR}/${MM10_TRANSCRIPTOME_FASTA_FILE}" > "${MURINE_TRANSCRIPTOME_DIR}/${MM10_TRANSCRIPTOME_FASTA}"
  rm "${MURINE_TRANSCRIPTOME_DIR}/${MM10_TRANSCRIPTOME_FASTA_FILE}"
  echo "Done Processing Mm10 Transcriptome"

  # Repeats Regions
  echo "Downloading MM10 B1 and B2 Repeats Table ${MM10_FTP_URL}${MM10_REGIONS_TABLE_FILE}"
  wget "${MM10_FTP_URL}${MM10_REGIONS_TABLE_FILE}"  --directory-prefix="${MURINE_REGIONS_DIR}"
  echo "Processing MM10  B1 and B2 Repeats Table ${MM10_REGIONS_TABLE_FILE}"
  zcat "${MURINE_REGIONS_DIR}/${MM10_REGIONS_TABLE_FILE}"| awk '{OFS ="\t"} (($13 ~/Alu/||$13 ~/^B2/) && $12 == "SINE"){print $6,$7,$8}' > "${MURINE_REGIONS_DIR}/${MM10_REGIONS_TABLE_FILE}.tmp"
  docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v "$(pwd)/${MURINE_REGIONS_DIR}":/data -t ${BEDTOOLS_DOCKER} /bin/bash -c "bedtools sort -i /data/${MM10_REGIONS_TABLE_FILE}.tmp | bedtools merge -i stdin | gzip > /data/${MM10_REGIONS_FILE}"
  rm "${MURINE_REGIONS_DIR}/${MM10_REGIONS_TABLE_FILE}" "${MURINE_REGIONS_DIR}/${MM10_REGIONS_TABLE_FILE}.tmp"
  echo "Done Processing MM10 B1 and B2 Repeats Table ${MM10_REGIONS_TABLE_FILE}"
  if [ "${DONT_GENERATE_GENOME_INDEXES}" = false ]
  then
    # JSD path is only genome name to match our Nextflow workflow requisits
    echo "Attempting to Create Genome Index of MM10 B1 and B2 Repeats ${MM10_REGIONS_FILE}"
    CMD="cd / ; /usr/bin/java -jar /bin/AEI/RNAEditingIndexer/lib/EditingIndexJavaUtils.jar GenerateIndex -i ${MURINE_REGIONS_DIR}/${MM10_REGIONS_FILE} -g ${MURINE_GENOME_DIR}/${MM10_GENOME_FASTA} -o ${MURINE_REGIONS_DIR}/${MM10_GENOME_FASTA}.${MM10_REGIONS_FILE}.GenomeIndex.jsd -b bedtools"
    # CMD="pwd"
    docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v "$(pwd)/${RESOURCES_DIR}":"/${RESOURCES_DIR}" -t ${RNA_EDITING_INDEX_DOCKER} /bin/bash -c "${CMD}"
    echo "Done Creating Genome Index of MM10 B1 and B2 Repeats ${MURINE_REGIONS_DIR}/${MM10_GENOME_FASTA}.${MM10_REGIONS_FILE}.GenomeIndex.jsd"
  fi

  # SNPs
  echo "Downloading MM10 Common Genomic SNPs Table ${MM10_FTP_URL}${MM10_SNPS_TABLE_FILE}"
  wget "${MM10_FTP_URL}${MM10_SNPS_TABLE_FILE}"  --directory-prefix="${MURINE_SNPS_DIR}"
  echo "Processing MM10 Genomic SNPs Table ${MM10_SNPS_TABLE_FILE}"
  zcat "${MURINE_SNPS_DIR}/${MM10_SNPS_TABLE_FILE}" | awk '{OFS ="\t"}($11=="genomic") {print $2,$3,$4,$7,$9,$10,$16,$25}'| gzip > "${MURINE_SNPS_DIR}/${MM10_SNPS_FILE}"
  rm "${MURINE_SNPS_DIR}/${MM10_SNPS_TABLE_FILE}"
  echo "Done Processing MM10 Common Genomic SNPs Table ${MM10_SNPS_TABLE_FILE}"

  # RefSeq Annotations
  echo "Downloading MM10 RefSeq Curated Table ${MM10_FTP_URL}${MM10_REFSEQ_TABLE_FILE}"
  wget "${MM10_FTP_URL}${MM10_REFSEQ_TABLE_FILE}"  --directory-prefix="${MURINE_REFSEQ_DIR}"
  echo "Processing MM10 RefSeq Curated Table ${MM10_REFSEQ_TABLE_FILE}"
  zcat "${MURINE_REFSEQ_DIR}/${MM10_REFSEQ_TABLE_FILE}"| awk '{OFS ="\t"} {print $3,$5,$6,$2,$13,$4,$10,$11}' |gzip > "${MURINE_REFSEQ_DIR}/${MM10_REFSEQ_FILE}"
  zcat "${MURINE_REFSEQ_DIR}/${MM10_REFSEQ_FILE}" | awk 'BEGIN {FS=OFS="\t"; print "Transcript\tGene"} { gsub(/\.[0-9_]+/,"",$4) ; print $4,$5}' > "${MURINE_REFSEQ_DIR}/${MM10_REFSEQ_FILE_FINAL}"
  rm "${MURINE_REFSEQ_DIR}/${MM10_REFSEQ_TABLE_FILE}"
  echo "Done Processing MM10 RefSeq Curated Table ${MM10_REFSEQ_TABLE_FILE}"

  echo "Downloading Mm10 RefSeq GTF Annotation File ${MM10_FTP_GENES_URL}${MM10_REFSEQ_GTF_ZIP_FILE}"
  curl "${MM10_FTP_GENES_URL}${MM10_REFSEQ_GTF_ZIP_FILE}" -o "${MURINE_REFSEQ_DIR}/${MM10_REFSEQ_GTF_ZIP_FILE}"
  echo "Saving Mm10 RefSeq Table Under: ${MURINE_REFSEQ_DIR}/${MM10_REFSEQ_GTF_FILE}"
  gunzip "${MURINE_REFSEQ_DIR}/${MM10_REFSEQ_GTF_ZIP_FILE}"
  echo "Done Processing Mm10 RefSeq GTF ${MM10_REFSEQ_GTF_FILE}"

  # Genes Expression
  echo "Downloading Murine Genes Expression Table ${MM10_GENE_EXPRESSION_FTP}${MM10_GENES_EXPRESSION_TABLE_FILE}"
  wget "${MM10_GENE_EXPRESSION_FTP}${MM10_GENES_EXPRESSION_TABLE_FILE}"  --directory-prefix="${MURINE_GENES_EXPRESSION_DIR}"
  echo "Processing MM10Genes Expression File ${MM10_GENES_EXPRESSION_TABLE_FILE}"
  docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v "$(pwd)/${MURINE_GENES_EXPRESSION_DIR}":/data -t ${BIGBED_TO_BED_DOCKER} /bin/bash -c "bigBedToBed data/${MM10_GENES_EXPRESSION_TABLE_FILE} data/${MM10_GENES_EXPRESSION_TABLE_FILE}.tmp"
  awk 'BEGIN{ FS=OFS="\t"}{print $1, $2, $3, $7, $9, $6}' "${MURINE_GENES_EXPRESSION_DIR}/${MM10_GENES_EXPRESSION_TABLE_FILE}.tmp" | gzip > "${MURINE_GENES_EXPRESSION_DIR}/${MM10_GENES_EXPRESSION_FILE}" 
  rm "${MURINE_GENES_EXPRESSION_DIR}/${MM10_GENES_EXPRESSION_TABLE_FILE}" "${MURINE_GENES_EXPRESSION_DIR}/${MM10_GENES_EXPRESSION_TABLE_FILE}.tmp"
  echo "Done Processing MM10 Genes Expression Table ${MM10_GENES_EXPRESSION_TABLE_FILE}"
fi

#---------------------------------------------------------------------------
# Create INI File
#---------------------------------------------------------------------------
if [ "${DONT_WRITE}" = false ]
then
  {
  echo "[DEFAULT]"
  echo "ResourcesDir = ${RESOURCES_DIR}"
  echo "BEDToolsDocker = ${BEDTOOLS_DOCKER}"
  echo "BigBEDToBEDDocker = ${BIGBED_TO_BED_DOCKER}"
  echo "RNAEditingIndexDocker = ${RNA_EDITING_INDEX_DOCKER}"

  echo "[hg38]"
  echo "Genome =  gs://${GCP_BUCKET_NAME}/${HUMAN_GENOME_DIR}/${HG38_GENOME_FASTA}"
  echo "RERegions = gs://${GCP_BUCKET_NAME}/${HUMAN_REGIONS_DIR}/${HG38_REGIONS_FILE}"
  echo "SNPs = gs://${GCP_BUCKET_NAME}/${HUMAN_SNPS_DIR}/${HG38_SNPS_FILE}"
  echo "RefSeq = gs://${GCP_BUCKET_NAME}/${HUMAN_REFSEQ_DIR}/${HG38_REFSEQ_FILE}"
  echo "RefSeqAll = gs://${GCP_BUCKET_NAME}/${HUMAN_REFSEQ_DIR}/${HG38_ALL_REFSEQ_FILE}"
  echo "GenesExpression = gs://${GCP_BUCKET_NAME}/${HUMAN_GENES_EXPRESSION_DIR}/${HG38_GENES_EXPRESSION_FILE}"
  echo ""

  echo "[mm10]"
  echo "Genome =  gs://${GCP_BUCKET_NAME}/${MURINE_GENOME_DIR}/${MM10_GENOME_FASTA}"
  echo "RERegions = gs://${GCP_BUCKET_NAME}/${MURINE_REGIONS_DIR}/${MM10_REGIONS_FILE}"
  echo "CDSRegions = gs://${GCP_BUCKET_NAME}/${MURINE_REGIONS_DIR}/${MM10_CDS_REGIONS_FILE}"
  echo "CDSRegionsRefSeqAll = gs://${GCP_BUCKET_NAME}/${MURINE_REGIONS_DIR}/${MM10_CDS_REFSEQ_ALL_REGIONS_FILE}"
  echo "SNPs = gs://${GCP_BUCKET_NAME}/${MURINE_SNPS_DIR}/${MM10_SNPS_FILE}"
  echo "RefSeq = gs://${GCP_BUCKET_NAME}/${MURINE_REFSEQ_DIR}/${MM10_REFSEQ_FILE}"
  echo "RefSeqAll = gs://${GCP_BUCKET_NAME}/${MURINE_REFSEQ_DIR}/${MM10_ALL_REFSEQ_FILE}"
  echo "GenesExpression = gs://${GCP_BUCKET_NAME}/${MURINE_GENES_EXPRESSION_DIR}/${MM10_GENES_EXPRESSION_FILE}"
  echo ""
   } > "${DBS_PATHS_INI}"

else
  echo "Info: Writing Resources.ini file is Disabled!"
fi