// nextflow - enable DSL 2
nextflow.enable.dsl=2

// region CheckParams
if(!params.srrACC_list){
    println "please set srrACC_list \n exiting..."
    exit 1
}
if ( params.isPaired == 'unset' || !params.isStranded == "unset") {
    println "please set isPaired & isStranded \n exiting..."
    exit 1
}
if(!params.organism){
    println "please set organism profile \n exiting..."
    exit 1
}
// endregion




// print current run parameters
log.info """\
         FROM SRRS to RNA-CEI  N F   P I P E L I N E
         =========================================================
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
            nextflow -c rna_editing.config -bg run rna_editing.nf -bucket-dir <nexflow_bucket_workdir> -profile <SE,stranded,RL75,hg38> --run_title <RUN_TITLE> --srrACC_list <SRR_LIST> > log.out 2> log.err &
            """
            .stripIndent()

}

if (params.help) {
    helpMessage()
    exit 0
}
// endregion




process download {
    cache 'lenient'
    cpus 2
    memory { task.attempt > 1 ? '16 GB' : '8 GB' }
    disk { task.attempt > 1 ? '200 GB' : '150 GB' }
    tag "download $srr"
    input:
        val srr
        val dirs_finished
        path "ngc_file"
    output:
        tuple val(srr), path('data/*fastq*'), emit: fastqs
        env fsize, emit: fastqsSize
    script:
        if(params.NGC_file){
            NGC_flag=" --ngc ./ngc_file"
        } else {
            NGC_flag=""
        }
        """
        # info
        echo "step download on $srr"
        # configure SRAtoolkit
        vdb-config --report-cloud-identity yes

        # fetch SRA file
        echo downloading
        # will not exit in case of error
        # because we will try download again
        {
            prefetch${NGC_flag} -C yes $srr -S no --max-size 100g &> download.log.txt 
        } || {
            echo "error downloading $srr"
            cat download.log.txt
        }
        echo "finished download first time (via cloud) $srr"
        # info
        # if *.sra file was downloaded via cloud then: extract FASTQ
        if [ -f ./${srr}/${srr}.sra ]; then
            echo "sra file was downloaded via cloud, extracting $srr ..."        
        else
            # remove bad *.sralite file if exist
            if [ -d ./${srr} ]; then
                echo "sra lite was downloaded, removing  $srr and downloading second time (not via cloud)"
                rm -r ./$srr
            fi
            echo "downloading second time (not via cloud) $srr"
            # fetch SRA file not via cloud
            vdb-config --report-cloud-identity no
            {
                prefetch${NGC_flag} -C yes $srr -S no --max-size 100g &> download.log.txt
            } || {
                echo "error downloading $srr"
                cat download.log.txt
            }
            echo "ls ."
            ls
            echo "ls dir"
            ls ${srr}
            # if file called something like SRR1234567_dbgap_#####.sra - rename it
            if compgen -G "./${srr}/${srr}*.sra" > /dev/null; then
                echo "files in dir:"
                ls ./${srr}
                # if ./${srr}/${srr}.sra not exist - rename the wrong file to the right name
                if [ ! -f ./${srr}/${srr}.sra ]; then
                    echo "renaming file"
                    mv ./${srr}/${srr}*.sra ./${srr}/${srr}.sra
                fi
            elif compgen -G "./${srr}*.sra" > /dev/null; then
                mkdir ${srr} || echo "dir exist"
                if [ ! -f ./${srr}/${srr}.sra ]; then
                    echo "renaming dbgap file"
                    mv ./${srr}*.sra ./${srr}/${srr}.sra
                fi
            fi
            echo "finished downloading second time (not via cloud) $srr"

        fi
        # extract FASTQ
        echo "extracting $srr"
        {
            fasterq-dump${NGC_flag} -e 2 -O ./data/ ./${srr}/${srr}.sra &> extract.log.txt

        } || {
            echo "failed extracting" 
            ls ./${srr}
            cat extract.log.txt
            echo "failed extracting" 
            exit 1
        }
        echo "finished extracting $srr"

        # delete unpaired file if exist:
        # if there is both _1  _2 files  delete third file
        if compgen -G "./data/${srr}*_1*" > /dev/null && compgen -G "./data/${srr}*_2*" > /dev/null; then
            if [ $params.isPaired = false ]; then
                echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! runing as SE but its PE $srr !!!!!!!!!!!!!!!!!!!!!"
            fi
            # if there is also unpaired file - remove it
            if [ -f ./data/${srr}.fastq ]; then
                # notify
                echo "removing third file"
                # remove
                rm ./data/${srr}.fastq
            fi
        else
            if [ $params.isPaired = true ]; then
                if [ -f ./data/${srr}.fastq ]; then
                    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! runinng as PE but its SE $srr !!!!!!!!!!!!!!!!!!!!!"
                else
                    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! no fastq was downloaded $srr !!!!!!!!!!!!!!!!!!!!!"
                fi
            fi

        fi

        echo -e "final downloaded files:\n\$(du -h ./data/*.fastq)"

        # send fastqs size for allocate disk for next steps:
        fsize=\$(du -B 1G ./data | tail -1 | cut -f1)
        # if fsize > params.fastqs_compression_threshold (GB) - compress fastq for being sure it will fit in fargate
        if (( \$fsize > ${params.fastqs_compression_threshold} )); then
            echo "compressing fastq files for $srr"
            gzip -f ./data/*.fastq
            # update fsize to compressed size
            gzipfsize=\$(du -B 1G ./data | tail -1 | cut -f1)
            echo "compressed fastq size is \$gzipfsize GB"
        else
            echo "fastq files are small enough, no need to compress"
        fi

        # notify
        echo "step download for $srr finished successfuly"
        """
}


process preprocess {
    cache 'lenient'
    disk { ds }
    memory '8 GB'
    cpus 4
    tag "fastp $srr"
    publishDir "${params.fastp_jsons}", mode: 'copy', overwrite: true, pattern: "*.json"
    input:
        tuple val(srr), path(fastqs)
        val fastqsSize
    output:
        tuple val(srr), path("outputs/*"), emit: fastqs
        path "${srr}.json", emit: json
    script:
        // dynamically compute disk size
        ds = fastqsSize as Integer
        // we will have 2 copies of the fastqs +10gb for resources and OS and docker
        ds = ((ds*3) + 10) as Integer
        // check it compatible with min&max of fargate
        if(ds < 21){
            ds=21
        }
        if (ds>200){
            println "$srr fastp data size is too big > 200GB"
            ds=200
        }
        ds="$ds GB"
        def read1 = params.isPaired ? "${fastqs[0]}" : "${fastqs}"
        def read2 = params.isPaired ? "${fastqs[1]}" : "NA"
        def fastp_read1_inputOutput = " -i ${read1} -o outputs/${read1}"
        def fastp_read2_inputOutput = params.isPaired ? " -I ${read2} -O outputs/${read2}" : ""
        def fastp_inputOutput = "${fastp_read1_inputOutput}${fastp_read2_inputOutput}"
        """
        # info
        echo "step preprocess on $srr"
        # create local directory
        mkdir outputs

        # info
        # don't trim polyG of reads -G
        # don't trim adapters -A
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
        fastp ${fastp_inputOutput} -j ${srr}.json -n 5 -q 25 -u 20 -e 30 --dont_eval_duplication -G -A --max_len1 \$(($params.fastp_read_length)) --length_required \$(($params.fastp_read_length-3)) -w 10

        
        echo -e "final processed files:\n\$(du -h outputs/*.fastq*)"

        # if no reads passed filters - files are empty
        if [ ! -s outputs/${read1} ] || ([ -f outputs/${read2} ] && [ ! -s outputs/${read2} ]); then
            # notify
            echo "no reads passed preprocessing for $srr"
            # stop process
            exit 1
        fi

        # notify
        echo "step preprocess for $srr finished successfuly"
        """
}

//quant genes expressions
process exp_quant {
    cache 'lenient'
    disk { ds }
    cpus 4
    memory '20 GB'
    tag "salmon $srr"
    publishDir "${params.genes_expressions}", mode: 'copy', overwrite: true, pattern: "*.sf"
    input:
        tuple val(srr), path(reads)
        val fastqsSize
        path transcripts_index
        path tx2id_geneMap 
    output:
        val finished, emit: finished
        path "*.sf", emit: SFs
    script:    
        finished = true
        // threads
        runThreadN = '8'
        // library type
        libType = 'A'
        // create command compatible for paired or single fastq
        def reads_commmand = params.isPaired ? "-1 ${reads[0]} -2 ${reads[1]}" : "-r ${reads}"
        // dynamically compute disk size 
        ds = fastqsSize as Integer
        // we need 40GB for salmon resources and OS and docker
        // also add some extra for temp files
        ds = (40 + ds*2) as Integer
        // check it compatible with min&max of fargate 
        if(ds < 21){
            ds=21
        }
        if (ds>200){
            println "$srr salmon data size is too big > 200GB"
            ds=200
        }
        ds="$ds GB"
        """
        # info
        echo "step exp_quant on $srr"
        # create local output directory
        salmon_dir=\${PWD}/salmon_genes_expression
        mkdir \$salmon_dir

        # log command
        # run salmon
        salmon quant -l ${libType} ${reads_commmand} -o \$salmon_dir -i $transcripts_index -g $tx2id_geneMap -p ${runThreadN}

        # rename output files:
        mv salmon_genes_expression/quant.sf ${srr}.quant.sf 
        mv salmon_genes_expression/quant.genes.sf ${srr}.quant.genes.sf
        
        # show output files
        ls salmon_genes_expression/

        # notify
        echo "step exp_quant for $srr finished successfuly"
        """
}


process alignment {
    cache 'lenient'
    disk { ds }
    cpus 8
    memory '40 GB'
    tag "star $srr"
    publishDir "${params.star_stats}", mode: 'copy', overwrite: true, pattern: "*.log"
    input:
        tuple val(srr), path(reads)
        val fastqsSize
        path genome_index_full
        val prev_finished
    output:
        tuple val(srr), path('outputs/*.bam'), emit: bam
        env bamS, emit: bamSize
        // will contain stats about the mapping
        path "${srr}_STAR_Final.log", emit: log
    script:
        // threads
        runThreadN = '32'
        // dynamically compute disk size
        fastqsSize = fastqsSize as Integer
        ds = fastqsSize
        // we need 40GB for STAR resources and OS and docker
        // it also wiil create bams with the size of the fastqs 
        ds = (50 + ds*2) as Integer
        // check it compatible with min&max of fargate 
        if(ds < 21){
            ds=21
        }
        if (ds>200){
            println "$srr star data size is too big > 200GB"
            ds=200
        }
        ds="$ds GB"
        // if its fastqs size is > 80 GB - it's should be compressed
        def read_command = fastqsSize > params.fastqs_compression_threshold ? "zcat" : "cat"
        // trying 'less' it should be able to handle both gzipped and unzipped files
        // def read_command = "less"
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
        STAR --readFilesCommand ${read_command} --readFilesIn ${reads} --genomeDir $genome_index_full --outSAMattributes All --outSAMtype BAM Unsorted --outFileNamePrefix outputs/$srr --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 600000 --outFilterMismatchNoverLmax 0.3 --outFilterMismatchNoverReadLmax 1 --outFilterMatchNminOverLread ${params.outFilterMatchNminOverLread} --outFilterMultimapNmax 1 --genomeLoad NoSharedMemory --runThreadN ${runThreadN}

        # info
        echo -e "BAM files:\n\$(du -h outputs/*.bam)" 
        # save bam size for output
        # in this docker -B isnt supported so we get size in mb
        bamS=\$(du -m outputs/*bam | tail -1 | cut -f1)
        # convert to gb and add 1 to ceil
        bamS=\$((\$bamS/1024+1))
        # rename and move stats to wd
        mv outputs/${srr}Log.final.out ${srr}_STAR_Final.log

        # notify
        echo "step alignment for $srr finished successfuly"
        """ 
}

process cmpileup {
    cache 'lenient'
    disk { ds }
    cpus '2'
    memory '16 GB'
    tag "cmpileup $srr"
    publishDir "${params.cmpileups_tables}", mode: 'copy', overwrite: true, pattern: "*.cmpileup"
    input:
        tuple val(srr), path(bam_file)
        val bamSize
        path genome_cmpileup_sites_index
        path cmpileup_ini
        path cmpileup_sites
        path genome_file
    output:
        val finished, emit: finished
        path "*.cmpileup", emit: cmpileup_res
        // path "*log", emit: cmp_logs
    script:
        // filter reads only if its paired end
        def filtering_commmand = params.isPaired ? "pe_filter=\\'-f\\ 2\\'" : ""
        finished = true
        // dynamically compute disk size
        ds = bamSize as Integer
        // the genome is 3G, the other file are tiny
        // add more 5 for OS and docker - give twice for being sure
        ds = (15 + ds*2) as Integer
        // check it compatible with min&max of fargate 
        if(ds < 21){
            ds=21
        }
        if (ds>200){
            println "$srr cmpileup data size is too big > 200GB"
            ds=200
        }
        ds="$ds GB"
        """
        echo "step cmpileup on $srr"
        # create dirs for results and logs
        mkdir results
        mkdir logs
        
        # save my currnet working directory
        Mywd=\$(pwd)
        # move bam file and genome index to CMP WD
        mv $bam_file /data/${srr}Aligned.sortedByCoord.out.bam
        echo "moving genome index"
        mkdir -p /data/output/cmpileups/
        mkdir -p /data/resources/
        mv $genome_cmpileup_sites_index /data/resources/$params.genome_cmpileup_sites_index_name
        mkdir -p /data/resources
        mv ${cmpileup_ini} /data/resources/${cmpileup_ini}
        mv ${cmpileup_sites} /data/resources/${cmpileup_sites}
        mv ${genome_file} /data/resources/${genome_file}
        cd /

        # run cmpileup
        python /home/biodocker/GGPS/Session/PipelineManger.py -t 0,7 -c /data/resources/${cmpileup_ini} -o /media -l /media/Log -d /data -f Aligned.sortedByCoord.out.bam --follow_links -a regions_coordinates_bed=\\'/data/resources/${cmpileup_sites}\\' genome_fasta=\\'/data/resources/${genome_file}\\' bam_file_suffix=\\'Aligned.sortedByCoord.out.bam\\' genome_index_path=\\'/data/resources/${params.genome_cmpileup_sites_index_name}\\' $filtering_commmand


        cd \$Mywd
        mv /media/${srr}/*cmpileup  .
        mv /media/Log/*log  .

        # notify
        echo "step cmpileup for $srr finished successfuly"
        """
}



process AEI_and_CEI {
    cache 'lenient'
    disk { ds }
    cpus 4
    memory { task.attempt >= 2 ? '22 GB' : '30 GB' }
    publishDir "${params.EI_results_path}", mode: 'copy', overwrite: true, pattern: "*.csv"
    publishDir "${params.EIrun_stats_path}", mode: 'copy', overwrite: true, pattern: "*.stats.txt.gz"
    tag "AEI and CEI $srr"
    input:
        // general input files
        tuple val(srr), path(bam_file)
        val bamSize
        // refseq annotation
        path refseq_file
        path expression_file
        path genome_file
        path snp_file
        // AEI input files
        path genome_regions_index_AEI
        path regions_AEI
        // CEI input files
        path genome_regions_index_CEI
        path regions_CEI
        val prev_finished
    output:
        path "*.csv", emit: results
        path "*.stats.txt", emit: stats
    script:
        // set EI according to stranded and PE parameters
        if(params.isPaired && params.isStranded)
            EI_stranded_paired_flag = ' --paired_end --stranded'
        else
            EI_stranded_paired_flag = ''
        def regions_index_path_override_alu = "-a genome_index_path=\\\'/data/sites_resources/${params.genome_regions_index_AEI_name}\\\'"
        def regions_index_path_override_CEI = "-a genome_index_path=\\\'/data/sites_resources/${params.genome_regions_index_CEI_name}\\\'"
        // dynamically compute disk size
        ds = bamSize as Integer
        // give enough space for the input,middle and results files
        ds = (21 + (ds*4)) as Integer
        ds=199
        // check it compatible with min&max of fargate 
        if(ds < 21){
            ds=21
        }
        if (ds>200){
            println "$srr aluAndCEI data size is too big > 200GB"
            ds=200
        }
        ds="$ds GB"
        """
        echo "step aluindex on $srr"
        # save path to stats file
        stat_file=\${PWD}/${srr}_runEI.stats.txt

        # create local directory
        mkdir -p /data/input/bams
        # create specific input paths to match RNAEditingIndexer tools
        mv $bam_file /data/input/bams/${srr}Aligned.sortedByCoord.out.bam
        echo "move region genome index"
        mkdir -p /data/output/cmpileups/
        mkdir -p /data/output/summary_dir
        mkdir -p /data/resources/
        mkdir -p /data/sites_resources/
        mv $genome_regions_index_AEI /data/sites_resources/${params.genome_regions_index_AEI_name}
        mv $regions_AEI /data/sites_resources/${regions_AEI}
        mv $refseq_file /data/resources/${refseq_file}
        mv $expression_file /data/resources/${expression_file}
        mv $genome_file /data/resources/${genome_file}
        mv $snp_file /data/resources/${snp_file}
        # save Nextflow workdir for future usage
        Mywd=\$(pwd)
        cd /

        echo "running AEI on $srr"
        # run AEI
        /usr/bin/time -v RNAEditingIndex -d /data/input/bams -l /data/output/logs_dir -o /data/output/cmpileups -os /data/output/summary_dir --follow_links --genome 'UserProvided' --refseq /data/resources/${refseq_file} --genes_expression /data/resources/${expression_file} -gf /data/resources/${genome_file} --snps /data/resources/${snp_file} -rb /data/sites_resources/${regions_AEI}${EI_stranded_paired_flag} $regions_index_path_override_alu 2> \$stat_file || cat \$stat_file

        echo "finished running AEI on $srr"

        # rename and move results
        cd \$Mywd
        # editing index summary
        mv /data/output/summary_dir/EditingIndex.csv ./${srr}_AluEditingIndex.csv
        echo "step aluindex for $srr finished successfuly"
        
        echo "step CEI on $srr"
        echo "removing AEI files"
        # remove all AEI files
        rm -r /data/output/logs_dir/* || echo "/data/output/logs_dir/ empty"
        rm -r /data/output/cmpileups/* || echo "/data/output/cmpileups/ empty"
        rm -r /data/output/summary_dir/* || echo "summary_dir empty"
        rm -r /data/sites_resources/* || echo "sites_resources empty"

        echo "move region and genome index"
        mv $genome_regions_index_CEI /data/sites_resources/${params.genome_regions_index_CEI_name}
        mv $regions_CEI /data/sites_resources/${regions_CEI}


        cd /
        # info 
        echo "running CEI on $srr"
        # run EI on CEI
        /usr/bin/time -v RNAEditingIndex -d /data/input/bams -l /data/output/logs_dir -o /data/output/cmpileups -os /data/output/summary_dir --follow_links --genome 'UserProvided' --refseq /data/resources/${refseq_file} --genes_expression /data/resources/${expression_file} -gf /data/resources/${genome_file} --snps /data/resources/${snp_file} -rb /data/sites_resources/${regions_CEI}${EI_stranded_paired_flag} $regions_index_path_override_CEI 2>> \$stat_file || cat \$stat_file

        echo "finished running CEI on $srr"

        # move and rename results
        cd \$Mywd
        mv /data/output/summary_dir/EditingIndex.csv ./${srr}_CEI.csv

        # notify
        echo "step CEI for $srr finished successfuly"
        """
}


// full without sorting bam
workflow {
    // convert srr list to channel
    srrs = Channel.fromPath(params.srrACC_list).splitText().map { it.replaceFirst(/\n/,'') }
    // download and proccess fastq files
    ngc_file = params.NGC_file ? params.NGC_file : "/dev/null"
    // download
    download(srrs,true,ngc_file)
    fastqs_size=download.out.fastqsSize
    // clean fastq files
    preprocess(download.out.fastqs,fastqs_size)
    cleaned_fastq_ch=preprocess.out.fastqs
    // quant gene expression
    if(params.runSalmon){
        exp_quant(cleaned_fastq_ch,fastqs_size,params.transcripts_index,params.tx2id_geneMap)
        expQ_finished = exp_quant.out.finished
    } else {
        expQ_finished = "not_running"
    }
    // create and proccess bams
    alignment(cleaned_fastq_ch,fastqs_size,params.genome_index_full,expQ_finished)
    bams_ch = alignment.out.bam
    bam_size = alignment.out.bamSize
    // if its human - run cmpileup on known editing CDS sites, and view results path
    cmpileup(bams_ch,bam_size,params.genome_cmpileup_sites_index,params.cmpileup_ini,params.cmpileup_sites,params.genome_file)
    CMP_finished = cmpileup.out.finished
    // Analize editing-index
    AEI_and_CEI(bams_ch,bam_size,params.refseq_file,params.expression_file,params.genome_file,params.snp_file,params.genome_regions_index_AEI,params.regions_AEI,params.genome_regions_index_CEI,params.regions_CEI,CMP_finished)
    AEI_and_CEI.out.results.view()
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


