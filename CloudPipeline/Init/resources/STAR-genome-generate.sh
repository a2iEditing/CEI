#---------------------------------------------------------------------------
# Init 
#---------------------------------------------------------------------------
RESOURCES_DIR=${1:-"Resources"};
STAR_INDEX_DIR=${2:-"STAR2.7.10b"};
NUM_THREADS=${3:-10};
STAR_docker=${4:-"levanonlab/star:2.7.10b--h9ee0642_0"}

#---------------------------------------------------------------------------
# Constants - see above script
#---------------------------------------------------------------------------
# paths
HUMAN="HomoSapiens"
MURINE="MusMusculus"
GENOME_DIR="Genomes"
REFSEQ_DIR="RefSeqAnnotations"
HG38_GENOME_FASTA="ucscHg38Genome.fa"
MM10_GENOME_FASTA="ucscMm10Genome.fa"
HG38_REFSEQ_FILE="hg38.ncbiRefSeq.gtf"
MM10_REFSEQ_FILE="mm10.ncbiRefSeq.gtf"

#---------------------------------------------------------------------------
# Index Genomes
#---------------------------------------------------------------------------

for readlen in $(seq 50 25 150); do 
	echo "Indexing hg38 genome for STAR with sjdbOverhang=${readlen} ..."
	# hg38
	path_to_references="${STAR_INDEX_DIR}/${HUMAN}/hg38/sjdbOverhang${readlen}"
	# create output directory
	[ ! -d $path_to_references ] && mkdir -p $path_to_references
	# get real paths of inputs and output
	path_to_references="$(realpath ${STAR_INDEX_DIR}/${HUMAN}/hg38/sjdbOverhang${readlen})"
	in_path_fasta="$(realpath ${RESOURCES_DIR}/${GENOME_DIR}/${HUMAN}/${HG38_GENOME_FASTA})"
	in_path_gtf="$(realpath ${RESOURCES_DIR}/${REFSEQ_DIR}/${HUMAN}/${HG38_REFSEQ_FILE})"
	# index genome
	docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v $path_to_references:/data -v $in_path_fasta:/inputs/Genome.fa -v $in_path_gtf:/inputs/ncbiRefSeq.gtf -t $STAR_docker /bin/bash -c "STAR --runMode genomeGenerate --genomeDir data --genomeFastaFiles /inputs/Genome.fa --sjdbGTFfile /inputs/ncbiRefSeq.gtf --sjdbOverhang ${readlen} --runThreadN ${NUM_THREADS} --outFileNamePrefix 'data/'"
	
	# mm10
	echo "Indexing mm10 genome for STAR with sjdbOverhang=${readlen} ..."
	path_to_references="${STAR_INDEX_DIR}/${MURINE}/mm10/sjdbOverhang${readlen}"
	# create output directory
	[ ! -d $path_to_references ] && mkdir -p $path_to_references
	# get real paths of inputs and output
	path_to_references="$(realpath ${STAR_INDEX_DIR}/${MURINE}/mm10/sjdbOverhang${readlen})"
	in_path_fasta="$(realpath ${RESOURCES_DIR}/${GENOME_DIR}/${MURINE}/${MM10_GENOME_FASTA})"
	in_path_gtf="$(realpath ${RESOURCES_DIR}/${REFSEQ_DIR}/${MURINE}/${MM10_REFSEQ_FILE})"
	# index genome
	docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v $path_to_references:/data -v $in_path_fasta:/inputs/Genome.fa -v $in_path_gtf:/inputs/ncbiRefSeq.gtf -t $STAR_docker /bin/bash -c "STAR --runMode genomeGenerate --genomeDir data --genomeFastaFiles /inputs/Genome.fa --sjdbGTFfile /inputs/ncbiRefSeq.gtf --sjdbOverhang ${readlen} --runThreadN ${NUM_THREADS} --outFileNamePrefix 'data/'"
done
