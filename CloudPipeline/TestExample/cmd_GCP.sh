BUCKET_NAME=
PROJECT_NAME=
REGION=
PLATFORM=GCP
NUM_THREADS=24

# Create working directory
mkdir -p work_test
cd work_test

# Clone git directory
git clone https://github.com/roni-fultheim/CEI.git

# Run initialization (long, requires at least 64G RAM)
# This example uses 24 CPUs in parallel
nohup sh $(pwd)/CEI/CloudPipeline/Init/init_main.sh $PLATFORM $BUCKET_NAME $NUM_THREADS > init.out 2> init.err &

# Run Nextflow command
outdir=$(pwd)/nextflow_output
[ ! -d $outdir ] && mkdir -p $outdir
cd $outdir
#export PATH="$PATH:$HOME/.local/bin"
nextflow -c $(pwd)/../CEI/CloudPipeline/GCP/SRA_pipeline/rna_editing.config -bg run $(pwd)/../CEI/CloudPipeline/GCP/SRA_pipeline/rna_editing.nf -profile hg38,PE,RL75,unstranded --project_name $PROJECT_NAME --region $REGION --bucket_name $BUCKET_NAME --run_title test_run_hg38_SE_RL75_unstranded --srrACC_list $(pwd)/../CEI/CloudPipeline/Test/hg38.PE.RL75.unstranded.txt > $outdir/run.out.txt 2> $outdir/run.err.txt &
