git clone https://github.com/roni-fultheim/CEI.git ./CEI
cd CEI
nohup sh CloudPipeline/Init/init_main.sh "AWS" <BUCKET_NAME> "24" > init.out 2> init.err &

nohup nextflow -c CloudPipeline/AWS/SRA_pipeline/rna_editing.config -bg run CloudPipeline/AWS/SRA_pipeline/rna_editing.nf -bucket-dir <nexflow_bucket_workdir> -profile <SE,stranded,RL75,hg38> --run_title <RUN_TITLE> --srrACC_list <SRR_LIST> > log.out 2> log.err &