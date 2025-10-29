#!/bin/bash
RESOURCES_DIR=${1:-"Resources"};
HG38_REGIONS_BED=${2:-"$(dirname $(readlink -f "$0"))/hg38.Alu3pUTR_minLen200_17022021.RepeatsInOppositeOrientationAtRegions.sorted.merged.bed.gz"};
MM10_REGIONS_BED=${3:-"$(dirname $(readlink -f "$0"))/mm10.AluB1AndB2_3pUTR_minLen100_17022021.RepeatsInOppositeOrientationAtRegions.sorted.merged.bed.gz"};
RNA_EDITING_INDEX_DOCKER='levanonlab/rna-editing-index-lite:1.0'
DONT_GENERATE_GENOME_INDEXES=${5:-false}

#---------------------------------------------------------------------------
# Constants - see above script
#---------------------------------------------------------------------------
HUMAN="HomoSapiens"
MURINE="MusMusculus"
GENOME_DIR="Genomes"
REGIONS_DIR="Regions"

HUMAN_GENOME_DIR="${RESOURCES_DIR}/${GENOME_DIR}/${HUMAN}"
HG38_GENOME_FASTA="ucscHg38Genome.fa"
MURINE_GENOME_DIR="${RESOURCES_DIR}/${GENOME_DIR}/${MURINE}"
MM10_GENOME_FASTA="ucscMm10Genome.fa"

HUMAN_REGIONS_DIR="${RESOURCES_DIR}/${REGIONS_DIR}/${HUMAN}"
MURINE_REGIONS_DIR="${RESOURCES_DIR}/${REGIONS_DIR}/${MURINE}"


#---------------------------------------------------------------------------
# HG38
#---------------------------------------------------------------------------
# create wanted filename
HG38_REGIONS_FILE="$(basename ${HG38_REGIONS_BED/RepeatsInOppositeOrientationAtRegions/InvertedRepeatsIn3pUTR})"
# copy to resource directory
cp ${HG38_REGIONS_BED} "${HUMAN_REGIONS_DIR}/${HG38_REGIONS_FILE}"

# Generate Indexes
if [ "${DONT_GENERATE_GENOME_INDEXES}" = false ]
then
    # JSD path is only genome name to match our Nextflow workflow requisits
    echo "Attempting to Create Genome Index of Hg38 Alu Inverted Repeats in 3'UTR ${HG38_REGIONS_FILE}"
    CMD="cd / ; /usr/bin/java -jar /bin/AEI/RNAEditingIndexer/lib/EditingIndexJavaUtils.jar GenerateIndex -i ${HUMAN_REGIONS_DIR}/${HG38_REGIONS_FILE} -g ${HUMAN_GENOME_DIR}/${HG38_GENOME_FASTA} -o ${HUMAN_REGIONS_DIR}/${HG38_GENOME_FASTA}.${HG38_REGIONS_FILE}.GenomeIndex.jsd -b bedtools"
    # CMD="pwd"
    docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v "./${RESOURCES_DIR}":"/${RESOURCES_DIR}" -t ${RNA_EDITING_INDEX_DOCKER} /bin/bash -c "${CMD}"
    echo "Done Creating Genome Index of Hg38 Alu Repeats ${HUMAN_REGIONS_DIR}/${HG38_GENOME_FASTA}.${HG38_REGIONS_FILE}.GenomeIndex.jsd"
fi

#---------------------------------------------------------------------------
# MM10
#---------------------------------------------------------------------------
# create wanted path
MM10_REGIONS_FILE="$(basename ${MM10_REGIONS_BED/RepeatsInOppositeOrientationAtRegions/InvertedRepeatsIn3pUTR})"
# copy to resource directory
cp ${MM10_REGIONS_BED} "${MURINE_REGIONS_DIR}/${MM10_REGIONS_FILE}"

if [ "${DONT_GENERATE_GENOME_INDEXES}" = false ]
then
    # JSD path is only genome name to match our Nextflow workflow requisits
    echo "Attempting to Create Genome Index of MM10 B1 and B2 Inverted Repeats in 3'UTR ${MM10_REGIONS_FILE}"
    CMD="cd / ; /usr/bin/java -jar /bin/AEI/RNAEditingIndexer/lib/EditingIndexJavaUtils.jar GenerateIndex -i ${MURINE_REGIONS_DIR}/${MM10_REGIONS_FILE} -g ${MURINE_GENOME_DIR}/${MM10_GENOME_FASTA} -o ${MURINE_REGIONS_DIR}/${MM10_GENOME_FASTA}.${MM10_REGIONS_FILE}.GenomeIndex.jsd -b bedtools"
    # CMD="pwd"
    docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v "./${RESOURCES_DIR}":"/${RESOURCES_DIR}" -t ${RNA_EDITING_INDEX_DOCKER} /bin/bash -c "${CMD}"
    echo "Done Creating Genome Index of MM10 B1 and B2 Repeats ${MURINE_REGIONS_DIR}/${MM10_GENOME_FASTA}.${MM10_REGIONS_FILE}.GenomeIndex.jsd"
fi