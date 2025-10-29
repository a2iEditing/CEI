// nextflow - enable DSL 2
nextflow.enable.dsl=2

// settings and configurations
params.srrACC_list = ''
params.bucket_adress=''
params.isStranded = ''
params.isPaired = ''

// delete BAM, FASTQ when no longer required to save storage
params.delete_big_files=true

// region CheckParams
if(!params.EI_results_path_relative || !params.fastp_jsons_relative){
    println "please set EI_results_path_relative & fastp_jsons_relative \n exiting..."
    exit 1
}

if(!params.bucket_adress  || !params.project_name || !params.mount_bucket_dir || !params.srrACC_list){
    println "please set bucket_adress & srrACC_list & mount_bucket_dir & project_name \n exiting..."
    exit 1
}


if ( params.isPaired == 'unset' || !params.isStranded == "unset") {
    println "please set isPaired & isStranded \n exiting..."
    exit 1
}

//check index params exist
if(!params.regions_AEI_relative || !params.genome_regions_index_AEI_name || !params.genome_regions_index_AEI_relative || !params.refseq_file_relative || !params.expression_file_relative || !params.genome_file_relative || !params.snp_file_relative ){
    println "please set Alu index parameters \n exiting..."
    exit 1
}

//check 3UTR params exist
if(!params.regions_3UTR_relative || !params.genome_regions_index_3UTR_name || !params.genome_regions_index_3UTR_relative){
    println "please set 3'UTR index parameters \n exiting..."
    exit 1
}

//check Salmon params exist
if(!params.transcripts_index_relative || !params.tx2id_geneMap_relative){
    println "please set 3'UTR index parameters \n exiting..."
    exit 1
}

// check star params exist
if(!params.genome_relative_path){
    println "please set genome_relative_path \n exiting..."
    exit 1
}
// endregion

// set AEI according to stranded and PE parameters
if(params.isPaired && params.isStranded)
    params.aei_stranded_paired_flag = '--paired_end --stranded'
else
    params.aei_stranded_paired_flag = ''


// print current run parameters
log.info """\
         FROM SRRS to RNA-CEI N F  P I P E L I N E
         =======================================================
         profile           : ${workflow.profile}
         result in         : ${params.EI_results_path}

        This project is a basic A-to-I RNA editing analysis workflow, designed for genome-wide largescale RNA editing analysis in human and mouse.
        The workflow can be easily adapted for other uses and organisms.

         """
         .stripIndent()


// region HelpMessage
def helpMessage() {
    log.info """\
            This project is a basic A-to-I RNA editing analysis workflow, designed for genome-wide largescale RNA editing analysis in human and mouse. The workflow can be easily adapted for other uses and organisms.

            FROM SRRS to RNA-CEI - HELP!!!!
            ===================================
            nextflow -c rna_editing.config -bg run rna_editing.nf -profile hg38,PE,RL75,unstranded --project_name $PROJECT_NAME --region $REGION --bucket_name $BUCKET_NAME --run_title test_run_hg38_SE_RL75_unstranded --srrACC_list CEI/CloudPipeline/Test/hg38.PE.RL75.unstranded.txt > $outdir/run.out.txt 2> $outdir/run.err.txt &
            """
            .stripIndent()

}

if (params.help) {
    helpMessage()
    exit 0
}
// endregion

//TODO in config!!!
/*
check for the withName A|B in process section - make docker for index look better

change CMP bed to orshaiUhillel.bed - ADDED, TODO after testing index and salmon

think if we want to separate the results directory of AEI and 3UTREI

*/


// create directories for the results 
process create_dirs {
    label 'stable'
    // move to label - next 4 lines
    machineType 'e2-small'
    cpus 2
    memory '2 GB'
    disk '10 GB'
    input:
        path proj_dir
    output:
        val finished
    script:
        finished = true
        """
        mkdir -p ${params.fastp_jsons}
        mkdir -p ${params.star_stats}
        mkdir -p ${params.EI_results_path}
        mkdir -p ${params.cmpileups_tables}
        mkdir -p ${params.genes_expressions}
        mkdir -p ${params.logs_path}
        """
}


process download {
    label 'error_prone'
    // move to label - next 4 lines
    // machineType 'e2-small'
    // cpus 2
    // memory '2 GB'


    machineType 'e2-standard-4'
    cpus 4
    memory '16 GB'
    disk '150 GB'

    tag "download $srr"
    input:
        path proj_dir
        val srr
        val dirs_finished
        path "ngc_file"
    output:
        tuple env(ds), val(srr), path('data/*fastq')
        // TODO check download succeed,remove unpaired file,check SE compatible
    script:
        if(params.NGC_file){
            NGC_flag=" --ngc ./ngc_file"
        } else {
            NGC_flag=""
        }
        def log_file = "${params.logs_path}/${srr}_runLog.txt"
        """
        # info
        echo "step download on $srr"
        echo "step download on $srr" > $log_file

        # configure SRAtoolkit
        vdb-config --report-cloud-identity yes

        # fetch SRA file
        echo downloading
        # will not exit in case of error
        # because we will try download again
        prefetch${NGC_flag} -C yes $srr -S no --max-size 100g || echo "error downloading $srr"
        echo "finished download first time $srr"
        echo "finished download first time (via cloud) $srr"

        # info
        # if *.sra file was downloaded via cloud then: extract FASTQ
        if [ -f ./${srr}/${srr}.sra ]; then
            echo "sra file was downloaded via cloud, extracting $srr ..."
            #echo "sra file was downloaded via cloud, extracting $srr ..." 
        else 
            # if *.sralite file was downloaded via cloud or other problem: download *.sra not via cloud and extract FASTQ
            # remove bad *.sralite file if exist
            if [ -d ./${srr} ]; then
                echo "sra lite was downloaded, removing  $srr and downloading second time (not via cloud)"
                rm -r ./$srr
            fi
            # fetch SRA file not via cloud
            vdb-config --report-cloud-identity no
            echo "downloading second time (not via cloud) $srr"
            prefetch${NGC_flag} -C yes $srr -S no --max-size 100g || echo "error downloading $srr" 
            echo "finished downloading second time (not via cloud) $srr"
            # info
        fi
        # extract FASTQ
        echo "extracting $srr"
        fasterq-dump${NGC_flag} -e 4 -O ./data/ ./${srr}/${srr}.sra || (echo "failed extracting" && exit 1)
        echo "finished extracting $srr"

        # delete unpaired file if exist:
        # if there is both _1  _2 files  delete third file
        if compgen -G "./data/${srr}*_1*" > /dev/null && compgen -G "./data/${srr}*_2*" > /dev/null; then
            if [ $params.isPaired = false ]; then
                echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! runing as SE but its PE $srr !!!!!!!!!!!!!!!!!!!!!"
                echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! runing as SE but its PE $srr !!!!!!!!!!!!!!!!!!!!!" >> $log_file
            fi
            # if there is also unpaired file - remove it
            if [ -f ./data/${srr}.fastq ]; then
                # notify
                echo "removing third file"
                echo "removing third file" >> $log_file
                # remove
                rm ./data/${srr}.fastq
            fi
        else
            if [ $params.isPaired = true ]; then
                if [ -f ./data/${srr}.fastq ]; then
                    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! runinng as PE but its SE $srr !!!!!!!!!!!!!!!!!!!!!"
                    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! runinng as PE but its SE $srr !!!!!!!!!!!!!!!!!!!!!" >> $log_file
                else
                    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! no fastq was downloaded $srr !!!!!!!!!!!!!!!!!!!!!"
                    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! no fastq was downloaded $srr !!!!!!!!!!!!!!!!!!!!!" >> $log_file
                fi
            fi

        fi

        # dynamically allocate disk for next steps:
        # set disk size we need double sized for next proccess(because fastp will create new fastqs)- you will also get 15gb more then you ask
        fsize=\$(du -B 1G ./data | tail -1 | cut -f1)
        echo "files size is \$fsize"
        echo -e "final downloaded files:\n\$(du -h ./data/*.fastq)" 
        ds=\$((\$fsize*2))
        ds="\$ds GB"

        # notify
        echo "step download for $srr finished successfuly"
        """
}

process preprocess {
    label 'stable'

    // move to label - next 4 lines
    machineType 'e2-highcpu-8'
    memory '8 GB'
    cpus 8
    disk { ds }

    tag "fastp $srr"
    input:
        path proj_dir
        tuple val(ds), val(srr), path('input_?.fastq')
    output:
        tuple val(srr), val(ds), path("outputs/*")
    script:
        def log_file = "${params.logs_path}/${srr}_runLog.txt"
        // todo update all resources for next stages
        """
        # info
        echo "step preprocess on $srr"
        # create local directory
        mkdir outputs

        # if sample is paired then: run paired-end code
        if [ $params.isPaired = true ]; then
            # info
            
            # PE inputs
            # PE outputs
            # don't trim polyG of reads
            # don't trim adapters
            # don't evaluate deduplication rates
            # number of threads
            # quality filtering: N bases per read <= 5
            # quality filtering: average quality per read >= 30
            # quality filtering: % of low-quality bases per read <= 20
            # quality filtering: low quality per base <= 25
            # length filtering: reads len >= params.read_length-10
            # length trimming: trim reads with len > (params.read_length+5), to (params.read_length+5)
            
            # log command
            # run
            fastp -i input_1.fastq -I input_2.fastq -o outputs/output_1.fastq -O outputs/output_2.fastq -j ${srr}.json -n 5 -q 25 -u 20 -e 30 --dont_eval_duplication --max_len1 \$(($params.read_length+5)) --length_required \$(($params.read_length-10)) -w 32
        else
            # info
            
            # if sample is paired then: run single-end code
            # SE inputs
            # SE outputs
            # don't trim polyG of reads
            # don't trim adapters
            # don't evaluate deduplication rates
            # number of threads
            # quality filtering: N bases per read <= 5
            # quality filtering: average quality per read >= 30
            # quality filtering: % of low-quality bases per read <= 20
            # quality filtering: low quality per base <= 25
            # length filtering: reads len >= params.read_length-10
            # length trimming: trim reads with len > (params.read_length+5), to (params.read_length+5)

            # run
            fastp -i input_1.fastq -o outputs/output_1.fastq -j ${srr}.json -n 5 -q 25 -u 20 -e 30 --dont_eval_duplication --max_len1 \$(($params.read_length+5)) --length_required \$(($params.read_length-10)) -w 32
        fi

        
        # copy JSON to remote results storage
        cp ${srr}.json ${params.fastp_jsons}/${srr}.json
        
        echo -e "final processed files:\n\$(du -h outputs/*.fastq)"

        # if no reads passed filters - files are empty
        if [ ! -s outputs/output_1.fastq ] || ([ -f outputs/output_2.fastq ] && [ ! -s outputs/output_2.fastq ]); then  
            # notify
            echo "no reads passed preprocessing for $srr"
            # stop process
            exit 1
        fi

        # notify
        echo "step preprocess for $srr finished successfuly"
        
        # remove FASTQ
        echo "removing original FASTQ file for $srr"
        rm \$(readlink -f input_1.fastq)
        if [ $params.isPaired = true ]; then
            echo "removing original P2 for $srr"
            rm \$(readlink -f input_2.fastq)
        fi
        """
}

//quant genes expressions
process exp_quant {
    label 'stable'
    
    // move to label - next 4 lines
    machineType 'e2-standard-8'
    cpus 8
    memory '32 GB'
    // TODO - selecet disk size for this step  - make sure you pass disk-size to the next step correctly
    disk { ds }

    tag "salmon $srr"
    input:
        path proj_dir
        tuple val(srr), val(ds), path(reads)
    output:
        val finished
    script:
        def log_file = "${params.logs_path}/${srr}_runLog.txt"
        finished = true
        // threads
        runThreadN = '8'
        // library type
        libType = 'A'
        // create command compatible for paired or single fastq
        def reads_commmand = params.isPaired ? "-1 ${reads[0]} -2 ${reads[1]}" : "-r ${reads}"
        """
        # info
        echo "step exp_quant on $srr"
        # create local output directory
        salmon_dir=\${PWD}/salmon_genes_expression
        mkdir \$salmon_dir

        # log command
        # run salmon
        salmon quant -l ${libType} ${reads_commmand} -o \$salmon_dir -i $params.transcripts_index -g $params.tx2id_geneMap -p ${runThreadN}

        # copy results to remote results storage
        cp \$salmon_dir/quant.sf ${params.genes_expressions}/${srr}.quant.sf 
        cp \$salmon_dir/quant.genes.sf ${params.genes_expressions}/${srr}.quant.genes.sf
        
        # notify
        echo "step exp_quant for $srr finished successfuly"
        """
}


process alignment {
    label 'stable'
    
    // move to label - next 4 lines
    machineType 'e2-highmem-8'
    cpus 8
    memory '64 GB'
    disk { ds }

    tag "star $srr"
    input:
        path proj_dir
        tuple val(srr), val(ds), path(reads)
        val prev_finished
    output:
        tuple val(srr), val(ds), path('outputs/*.bam')
    script:
        def log_file = "${params.logs_path}/${srr}_runLog.txt"
        // threads
        runThreadN = '8'
        
        // dynamically compute disk size
        ds = ds.split(' ')[0] as Integer
        ds = (20 + ds*(3/5)) as Integer
        ds="$ds GB"
        """
        # info
        echo "step alignment on $srr"
        # create local directory
        mkdir outputs
        
        # log command
        # run star
        # input parameters
        # output parameters
        # algorithm parameters
        # running resources
        STAR --readFilesCommand cat --readFilesIn ${reads} --genomeDir ${params.genome_index_full} --outSAMattributes All --outSAMtype BAM Unsorted --outFileNamePrefix outputs/$srr --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 600000 --outFilterMismatchNoverLmax 0.3 --outFilterMismatchNoverReadLmax 1 --outFilterMatchNminOverLread  0.66 --outFilterMultimapNmax 1 --genomeLoad NoSharedMemory --runThreadN ${runThreadN}

        # info
        echo -e "BAM files:\n\$(du -h outputs/*.bam)" 

        # copy run statistics to remote results storage
        cp outputs/${srr}Log.final.out ${params.star_stats}/${srr}_STAR_Final.log
        # copy run log to remote log storage
        cp outputs/${srr}Log.out ${params.logs_path}/${srr}_STAR.log

        # notify
        echo "step alignment for $srr finished successfuly"

        # remove FASTQ
        echo "removing processed FASTQ file for $srr"
        rm \$(readlink -f ${reads[0]})
        if [ $params.isPaired = true ]; then
            echo "removing processed P2 for $srr"
            rm \$(readlink -f ${reads[1]})
        fi
        """ 
}

process cmpileup {
    label 'stable'

    // move to label - next 4 lines
    machineType 'e2-standard-2'
    disk { ds }
    cpus 2
    memory '8 GB'

    tag "cmpileup $srr"
    input:
        path proj_dir
        tuple val(srr), val(ds), path(bam_file)
    output:
        val path_res, emit: cmpileup_res
        // path "/media", emit: cmp_logs
    script:
        def log_file = "${params.logs_path}/${srr}_runLog.txt"
        path_res="${params.cmpileups_tables}/${srr}*cmpileup"
        // filter reads only if its paired end
        def filtering_commmand = params.isPaired ? "pe_filter=\\'-f\\ 2\\'" : ""
        """
        echo "step cmpileup on $srr"
        # create dirs for results and logs
        mkdir results
        mkdir logs
        
        # save my currnet working directory
        Mywd=\$(pwd)
        # move bam file to CMP WD
        mv $bam_file /data/${srr}Aligned.sortedByCoord.out.bam
        cd /

        echo "link to sites genome index"
        mkdir -p /data/output/cmpileups/
        ln -s $params.genome_cmpileup_sites_index /data/output/cmpileups/$params.genome_cmpileup_sites_index_name

        # log command
        # run cmpileup  
        python /home/biodocker/GGPS/Session/PipelineManger.py -t 0,7 -c ${params.cmpileup_ini} -o /media -l /media/Log -d /data -f Aligned.sortedByCoord.out.bam --follow_links -a regions_coordinates_bed=\\'${params.cmpileup_sites}\\' genome_fasta=\\'${params.genome_file}\\' bam_file_suffix=\\'Aligned.sortedByCoord.out.bam\\' genome_index_path=\\'/data/output/cmpileups/${params.genome_cmpileup_sites_index_name}\\' $filtering_commmand

        # copy results to result dir
        mkdir -p ${params.cmpileups_tables}

        cp /media/${srr}/*cmpileup  ${params.cmpileups_tables}


        
        # notify
        echo "step cmpileup for $srr finished successfuly"
        # echo -e "step cmpileup for $srr finished successfuly" >> $log_file
        """
}

process aluAnd3UTR_Eindex {
    label 'stable'

    // move to label - next 4 lines
    machineType 'e2-highmem-2'
    disk { ds }
    cpus 2
    memory '16 GB'

    tag "alu & 3UTR index $srr"
    input:
        path proj_dir
        tuple val(srr), val(ds), path(bam_file)
        val prev_finished
    output:
        val (path_res_alu)
        val (path_res_3UTR) 
    script:
        def log_file = "${params.logs_path}/${srr}_runLog.txt"
        path_res_alu="${params.EI_results_path}/${srr}_AluEditingIndex.csv"
        path_res_3UTR="${params.EI_results_path}/${srr}_3UTREditingIndex.csv"
        def regions_index_path_override_alu = "-a genome_index_path=\\\'/data/output/cmpileups/${params.genome_regions_index_AEI_name}\\\'"
        def regions_index_path_override_3UTR = "-a genome_index_path=\\\'/data/output/cmpileups/${params.genome_regions_index_3UTR_name}\\\'"        
                """
        # info
        echo "step aluindex on $srr"
        # create local directory
        mkdir -p /data/input/bams
        # create specific input paths to match RNAEditingIndexer tools
        mv $bam_file /data/input/bams/${srr}Aligned.sortedByCoord.out.bam
        # save Nextflow workdir for future usage
        Mywd=\$(pwd)
        cd /

        echo "link to region genome index"
        mkdir -p /data/output/cmpileups/
        ln -s $params.genome_regions_index_AEI /data/output/cmpileups/$params.genome_regions_index_AEI_name

        # info 
        echo "running AEI on $srr"
        # log command
        # run AEI
        RNAEditingIndex -d /data/input/bams -l /data/output/logs_dir -o /data/output/cmpileups -os \$Mywd/summary_dir --follow_links --genome 'UserProvided' --refseq $params.refseq_file --genes_expression $params.expression_file -gf $params.genome_file --snps $params.snp_file -rb $params.regions_AEI $params.aei_stranded_paired_flag $regions_index_path_override_alu
        
        echo "finished running AEI on $srr"

        # copy CSV to remote results storage
        cd \$Mywd
        cp \$Mywd/summary_dir/EditingIndex.csv ${params.EI_results_path}/${srr}_AluEditingIndex.csv
        # copy run log to remote log storage
        # cp /data/output/logs_dir/*log ${params.logs_path}/${srr}_AEI.log
        
        # notify
        echo "step aluindex for $srr finished successfuly"
        
        echo "step 3UTREditingIndex on $srr"
        cd /
        # remove all AEI files
        rm -r /data/output/logs_dir/* || echo "/data/output/logs_dir/ empty"
        rm -r /data/output/cmpileups/* || echo "/data/output/cmpileups/ empty"
        rm -r \$Mywd/summary_dir/* || echo "summary_dir empty"

          

        echo "link to region genome index"
        mkdir -p /data/output/cmpileups/
        ln -s $params.genome_regions_index_3UTR /data/output/cmpileups/$params.genome_regions_index_3UTR_name

        # info 
        echo "running 3UTREditingIndex on $srr"
        # run EI on 3UTR
        RNAEditingIndex -d /data/input/bams -l /data/output/logs_dir -o /data/output/cmpileups -os \$Mywd/summary_dir --follow_links --genome 'UserProvided' --refseq $params.refseq_file --genes_expression $params.expression_file -gf $params.genome_file --snps $params.snp_file -rb $params.regions_3UTR $params.aei_stranded_paired_flag $regions_index_path_override_3UTR
        
        echo "finished running 3UTREditingIndex on $srr"

        # copy CSV to remote results storage
        cd \$Mywd
        cp \$Mywd/summary_dir/EditingIndex.csv ${params.EI_results_path}/${srr}_3UTREditingIndex.csv
        # copy run log to remote log storage
        # cp /data/output/logs_dir/*log ${params.logs_path}/${srr}_EI3UTR.log

        # notify
        echo "step 3UTREditingIndex for $srr finished successfuly"

        # remove bams
        rm \$(readlink -f /data/input/bams/${srr}Aligned.sortedByCoord.out.bam)
        """
}



// full without sorting bam
workflow {
    // create results dirs
    create_dirs_finished = create_dirs(params.project_dir).collect() // collect() for having value channel
    // convert srr list to channel
    srrs = Channel.fromPath(params.srrACC_list).splitText().map { it.replaceFirst(/\n/,'') }
    // download and proccess fastq files
    ngc_file = params.NGC_file ? params.NGC_file : "/dev/null"
    downloaded_fastq = download(params.project_dir,srrs,create_dirs_finished,ngc_file)
    clean_fastq_ch = preprocess(params.project_dir,downloaded_fastq)
    // quant gene expression
    // note, clean_fastq_ch==clean_fastq_ch_afterQuant, but letting me run the processes serialy
    expQ_finished = exp_quant(params.project_dir,clean_fastq_ch)
    // create and proccess bams
    bams = alignment(params.project_dir,clean_fastq_ch, expQ_finished)
    // if its human - run cmpileup on known editing CDS sites, and view results path
    cmpileup(params.project_dir,bams)
    CMP_finished = cmpileup.out.cmpileup_res

    // Analyze editing-index
    (AEI_results, UTR3EI_finished) = aluAnd3UTR_Eindex(params.project_dir,bams,CMP_finished)
    UTR3EI_finished.view()
}

workflow.onComplete {
    log.info("""
    Complete: Workflow Ended
    """)
}

workflow.onError {
    log.info("""
    Error: Workflow Ended With Error
    """)
}


