#!/bin/bash
RESOURCES_DIR=${1:-"Resources"};
SITES_FILE_HUMAN=${2:-"$(dirname $(readlink -f "$0"))/CDS_Sites/gabay_etal_natureCom_AG_sites_pmid35246538.csv"};
SITES_FILE_MOUSE=${3:-"$(dirname $(readlink -f "$0"))/Conserved_Sites/pinto_etal_genomeBiol_conserved_sites_pmid24393560.liftOver_from_hg19_to_mm10.bed"};
RNA_EDITING_INDEX_DOCKER=${4:-'levanonlab/rna-editing-index-lite:1.0'}
DONT_GENERATE_GENOME_INDEXES=${5:-false}

#---------------------------------------------------------------------------
# Constants - see above script
#---------------------------------------------------------------------------
# paths
HUMAN="HomoSapiens"
MOUSE="MusMusculus"

GENOME_DIR="Genomes"
SITES_DIR="Sites"

HUMAN_GENOME_DIR="${RESOURCES_DIR}/${GENOME_DIR}/${HUMAN}"
HG38_GENOME_FASTA="ucscHg38Genome.fa"
MOUSE_GENOME_DIR="${RESOURCES_DIR}/${GENOME_DIR}/${MOUSE}"
MM10_GENOME_FASTA="ucscMm10Genome.fa"

HUMAN_SITES_DIR="${RESOURCES_DIR}/${SITES_DIR}/${HUMAN}"
HG38_SITES="hg38Sites_GabayEtAl_NatureCom_PMID35246538.AG.bed"
MOUSE_SITES_DIR="${RESOURCES_DIR}/${SITES_DIR}/${MOUSE}"
MM10_SITES="mm10Sites_PintoEtAl_GenomeBiol_PMID24393560.bed"

echo $SITES_FILE
mkdir -p "${HUMAN_SITES_DIR}"
mkdir -p "${MOUSE_SITES_DIR}"

# create human BED file
cat "${SCRIPT_DIR}/${SITES_FILE_HUMAN}" | awk 'BEGIN { FS = "," } ;NR>1 {OFS="\t"; print $1,$2-1,$2}' > "${HUMAN_SITES_DIR}/${HG38_SITES}"
# copy mouse BED file
cp ${SITES_FILE_MOUSE} ${MOUSE_SITES_DIR}/${MM10_SITES}

#---------------------------------------------------------------------------
# Generate Indexes
#---------------------------------------------------------------------------
if [ "${DONT_GENERATE_GENOME_INDEXES}" = false ]
then
    # JSD path is only genome name to match our Nextflow workflow prequisits
    echo "Attempting to Create Genome Index of Hg38 CMpileup Sites "${HUMAN_SITES_DIR}/${HG38_SITES}""
    CMD="cd / ; /usr/bin/java -jar /bin/AEI/RNAEditingIndexer/lib/EditingIndexJavaUtils.jar GenerateIndex -i ${HUMAN_SITES_DIR}/${HG38_SITES} -g ${HUMAN_GENOME_DIR}/${HG38_GENOME_FASTA} -o ${HUMAN_SITES_DIR}/${HG38_GENOME_FASTA}.${HG38_SITES}.GenomeIndex.jsd -b bedtools"
    # CMD="pwd"
    docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v "./${RESOURCES_DIR}":"/${RESOURCES_DIR}" -t ${RNA_EDITING_INDEX_DOCKER} /bin/bash -c "${CMD}"
    echo "Done Creating Genome Index of Hg38 CMpileup Sites ${HUMAN_SITES_DIR}/${HG38_GENOME_FASTA}.${HG38_SITES}.GenomeIndex.jsd"

    # JSD path is only genome name to match our Nextflow workflow prequisits
    echo "Attempting to Create Genome Index of Mm10 CMpileup Sites "${MOUSE_SITES_DIR}/${MM10_SITES}""
    CMD="cd / ; /usr/bin/java -jar /bin/AEI/RNAEditingIndexer/lib/EditingIndexJavaUtils.jar GenerateIndex -i ${MOUSE_SITES_DIR}/${MM10_SITES} -g ${MOUSE_GENOME_DIR}/${MM10_GENOME_FASTA} -o ${MOUSE_SITES_DIR}/${MM10_GENOME_FASTA}.${MM10_SITES}.GenomeIndex.jsd -b bedtools"
    # CMD="pwd"
    docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v "./${RESOURCES_DIR}":"/${RESOURCES_DIR}" -t ${RNA_EDITING_INDEX_DOCKER} /bin/bash -c "${CMD}"
    echo "Done Creating Genome Index of Mm10 CMpileup Sites ${MOUSE_SITES_DIR}/${MM10_GENOME_FASTA}.${MM10_SITES}.GenomeIndex.jsd"
fi