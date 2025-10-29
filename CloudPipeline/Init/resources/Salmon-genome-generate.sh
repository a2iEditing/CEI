#!/usr/bin/env bash

RESOURCES_DIR=${1:-"Resources"};
SALMON_INDEX_DIR=${2:-"Salmon1.10.2"};
NUM_THREADS=${3:-14};
SALMON_DOCKER=${4:-"levanonlab/salmon:1.10.2--hecfa306_0"}

#---------------------------------------------------------------------------
# Constants - see above script
#---------------------------------------------------------------------------
# paths
HUMAN="HomoSapiens"
MURINE="MusMusculus"
GENOME_DIR="Genomes"
HG38_GENOME_FASTA="ucscHg38Genome.fa"
MM10_GENOME_FASTA="ucscMm10Genome.fa"
TRANSCRIPTOME_DIR="Transcriptomes"
HG38_TRANSCIPTOME_FASTA="ucscHg38RefSeqTranscriptome.fa"
MM10_TRANSCIPTOME_FASTA="ucscMm10RefSeqTranscriptome.fa"

#---------------------------------------------------------------------------
# HG38
#---------------------------------------------------------------------------
path_to_references="${SALMON_INDEX_DIR}/${HUMAN}/hg38"
# create output directory
[ ! -d $path_to_references ] && mkdir -p $path_to_references
# get real paths of inputs and output
path_to_references="$(realpath ${SALMON_INDEX_DIR}/${HUMAN}/hg38)"
in_path_genome="$(realpath ${RESOURCES_DIR}/${GENOME_DIR}/${HUMAN}/${HG38_GENOME_FASTA})"
in_path_transcriptome="$(realpath ${RESOURCES_DIR}/${GENOME_DIR}/${HUMAN}/${HG38_TRANSCIPTOME_FASTA})"

DECOYS_FILENAME="${HG38_GENOME_FASTA}.decoys.txt"
GENTROME_FILENAME="${HG38_TRANSCIPTOME_FASTA}.${HG38_GENOME_FASTA}.gentrome.fa.gz"

# TODO delete
# print all commands
set -x

# Salmon indexing requires the names of the genome targets, which is extractable by using the grep command
grep "^>" ${RESOURCES_DIR}/${GENOME_DIR}/${HUMAN}/${HG38_GENOME_FASTA} | cut -d " " -f 1 | sed -e 's/>//g' > $path_to_references/${DECOYS_FILENAME}

# Along with the list of decoys salmon also needs the concatenated transcriptome and genome reference file for index.
# NOTE: the genome targets (decoys) should come after the transcriptome targets in the reference
cat ${RESOURCES_DIR}/${TRANSCRIPTOME_DIR}/${HUMAN}/${HG38_TRANSCIPTOME_FASTA} ${RESOURCES_DIR}/${GENOME_DIR}/${HUMAN}/${HG38_GENOME_FASTA} > $path_to_references/${GENTROME_FILENAME}

# index
docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v $path_to_references:/data -t $SALMON_DOCKER /bin/bash -c "salmon index -t /data/${GENTROME_FILENAME} -d /data/${DECOYS_FILENAME} -i /data -p ${NUM_THREADS}"

#---------------------------------------------------------------------------
# MM10
#---------------------------------------------------------------------------
# TODO delete
# don't print commands
set +x

path_to_references="${SALMON_INDEX_DIR}/${MURINE}/mm10"
# create output directory
[ ! -d $path_to_references ] && mkdir -p $path_to_references
# get real paths of inputs and output
path_to_references="$(realpath ${SALMON_INDEX_DIR}/${MURINE}/mm10)"
in_path_genome="$(realpath ${RESOURCES_DIR}/${GENOME_DIR}/${MURINE}/${MM10_GENOME_FASTA})"
in_path_transcriptome="$(realpath ${RESOURCES_DIR}/${GENOME_DIR}/${MURINE}/${MM10_TRANSCIPTOME_FASTA})"

DECOYS_FILENAME="${MM10_GENOME_FASTA}.decoys.txt"
GENTROME_FILENAME="${MM10_TRANSCIPTOME_FASTA}.${MM10_GENOME_FASTA}.gentrome.fa.gz"

# TODO delete
# print all commands
set -x

# Salmon indexing requires the names of the genome targets, which is extractable by using the grep command
grep "^>" ${RESOURCES_DIR}/${GENOME_DIR}/${MURINE}/${MM10_GENOME_FASTA} | cut -d " " -f 1 | sed -e 's/>//g' > $path_to_references/${DECOYS_FILENAME}

# Along with the list of decoys salmon also needs the concatenated transcriptome and genome reference file for index.
# NOTE: the genome targets (decoys) should come after the transcriptome targets in the reference
cat ${RESOURCES_DIR}/${TRANSCRIPTOME_DIR}/${MURINE}/${MM10_TRANSCIPTOME_FASTA} ${RESOURCES_DIR}/${GENOME_DIR}/${MURINE}/${MM10_GENOME_FASTA} > $path_to_references/${GENTROME_FILENAME}

# index
docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v $path_to_references:/data -t $SALMON_DOCKER /bin/bash -c "salmon index -t /data/${GENTROME_FILENAME} -d /data/${DECOYS_FILENAME} -i /data -p ${NUM_THREADS}"
