// --- params check:
if(!params.gdc_UUID_list){
    println "please set UUID_list \n exiting..."
    exit 1
}

params.isStranded = 'unset'
params.isPaired = 'unset'
if ( params.isPaired == 'unset' || !params.isStranded == "unset") {
    println "please set isPaired & isStranded \n exiting..."
    exit 1
}

//check index params exist
if(!params.regions_EI_relative || !params.genome_regions_index_EI_name || !params.genome_regions_index_EI_relative || !params.refseq_file_relative || !params.expression_file_relative || !params.genome_file_relative || !params.snp_file_relative ){
    println "please set Alu index parameters \n exiting..."
    exit 1
}

// set EI according to stranded and PE parameters
if(params.isPaired && params.isStranded)
    params.ei_stranded_paired_flag = '--paired_end --stranded'
else
    params.ei_stranded_paired_flag = ''



// print current run parameters
def index_type = params.CEI ? 'CytoplasmicEditingIndex' : 'AluEditingIndex'
log.info """\
         TCGA RNA-CEI N F   P I P E L I N E
         =========================================================
         Index type        : ${index_type}
         Profile           : ${workflow.profile}
         Result in         : ${params.EI_results_path}

         """
         .stripIndent()


// region HelpMessage
def helpMessage() {
    log.info """\
            This project is a basic A-to-I RNA editing analysis workflow, designed for genome-wide largescale RNA editing analysis in TCGA. The workflow can be easily adapted for other uses and organisms.
        
            TCGA RNA-CEI - HELP!!!!
            ===================================
            first you must set all user parameters in the TCGA_rna_editing.awsFargate.user_params config file.
            run the pipeline with the following command:
            run with:
            nextflow -bg -c TCGA_rna_editing.awsFargate.config run TCGA_rna_editing.awsFargate.nf --run_title <title> --gdc_UUID_list <UUIDs_list_file_path> --CEI <true/false> -profile <PE/SE>,<stranded/unstranded> > log.out 2> log.err &
            """
            .stripIndent()

}

if (params.help) {
    helpMessage()
    exit 0
}



// download bam files from GDC
process download {
    cache 'lenient'
    label 'run_batch'
    label 'process_low'

    tag "download $gdc_UUID"
    input:
        val gdc_UUID
        path GDC_token
        val dir_created
    output:
        tuple val(gdc_UUID), path("./${gdc_UUID}/*.bam"), emit: bam_file
        env fsize, emit: bamsize
        // path "${gdc_UUID}.gdc.out.txt"
    script:
        """
        # info
        echo "step download $gdc_UUID"
        gdc-client download $gdc_UUID -t $GDC_token &> gdc.out.txt || echo "download failed"
        echo "download finished"
        echo "bam size:"
        du -B 1G ./${gdc_UUID} || cat gdc.out.txt
        bam_count=\$(ls -1 ./${gdc_UUID}/*.bam 2>/dev/null | wc -l)
        if [ \$bam_count == 0 ]; then 
            echo ERROR_NO_BAM
            exit 1
        fi 
        mv gdc.out.txt ${gdc_UUID}.gdc.out.txt
        fsize=\$(du -B 1G ./${gdc_UUID}/*.bam | tail -1 | cut -f1)
        """

}

// run Alu index
process aluindex {
    cache 'lenient'
    label 'run_batch'
    label 'process_medium'
    disk { ds }
    publishDir "${params.EI_results_path}", mode: 'copy', overwrite: true, pattern: "*.csv"
    tag "alu index $gdc_UUID"
    input:
        // general input files
        tuple val(gdc_UUID), path(bam_file)
        val bamSize
        // refseq annotation
        path refseq_file
        path expression_file
        path genome_file
        path snp_file
        // EI input files
        path genome_regions_index_EI
        path regions_EI
    output:
        path "*.csv", emit: results
    script:
        def regions_index_path_override_alu = "-a genome_index_path=\\\'/data/sites_resources/${params.genome_regions_index_EI_name}\\\'"
        // dynamically compute disk size
        ds = bamSize as Integer
        // give twice the bam size for tmp files/cmpileups+filterd bam
        // refseq annotation is 1G, the genome is 3G and other files are small
        // also we want 2G for the cmpileup tables
        // add more 5G for OS and docker
        ds = (20 + (ds*2)) as Integer
        // check it compatible with min&max of fargate 
        if(ds < 21){
            ds=21
        }
        if (ds>200){
            println "data size is too big > 200GB"
            ds=200
        }
        ds="$ds GB"
        """
        # info
        echo "step aluindex on $gdc_UUID"
        # create local directory
        mkdir -p /data/input/bams

        ###### filter bam - !!! we also use it to transfer the Bam file !!!
        echo "filtering bam"
        samtools view -h -b -q 255 -F 2308 -f 2 -@ 10 -o /data/input/bams/${gdc_UUID}Aligned.sortedByCoord.out.bam $bam_file
        rm $bam_file
        echo "filtering done"
        #######

        # check if bam is not empty - probably SE
        num_reads=\$(samtools view -c /data/input/bams/${gdc_UUID}Aligned.sortedByCoord.out.bam)

        # if bam is not empty
        if (( \$num_reads >  0 )); then

            # create specific input paths to match RNAEditingIndexer tools
            echo "moving Eindex resources"
            # move other files
            mkdir -p /data/output/cmpileups/
            mkdir -p /data/output/summary_dir
            mkdir -p /data/resources/
            mkdir -p /data/sites_resources/
            mv $genome_regions_index_EI /data/sites_resources/${params.genome_regions_index_EI_name}
            mv $regions_EI /data/sites_resources/${regions_EI}
            mv $refseq_file /data/resources/${refseq_file}
            mv $expression_file /data/resources/${expression_file}
            mv $genome_file /data/resources/${genome_file}
            mv $snp_file /data/resources/${snp_file}
            # save Nextflow workdir for future usage
            Mywd=\$(pwd)
            cd /


            # info 
            echo "running EI on $gdc_UUID"
            # log command
            # run EI
            RNAEditingIndex -d /data/input/bams -l /data/output/logs_dir -o /data/output/cmpileups -os /data/output/summary_dir --follow_links --genome 'UserProvided' --refseq /data/resources/${refseq_file} --genes_expression /data/resources/${expression_file} -gf /data/resources/${genome_file} --snps /data/resources/${snp_file} -rb /data/sites_resources/${regions_EI} $params.ei_stranded_paired_flag ${regions_index_path_override_alu}
            
            echo "finished running EI on $gdc_UUID"

            # rename and move results
            cd \$Mywd
            mv /data/output/summary_dir/EditingIndex.csv ./${gdc_UUID}_AluEditingIndex.csv

        # if bam is empty - probably SE
        else
            # notify
            echo "bam is empty"   
            # create empty files with informative names for the next stages - we dont want the pipeline to crash
            touch ./${gdc_UUID}_AluEditingIndex.EMPTY.csv
            touch ./${gdc_UUID}_EI.EMPTY.log
        fi
        # notify
        echo "step aluindex for $gdc_UUID finished"
        """
}


// full without sorting bam
workflow {
    // create_dirs_finished =true
    gdc_UUIDs = Channel.fromPath(params.gdc_UUID_list, checkIfExists: true).splitText().map { it.replaceFirst(/\n/,'') }
    download(gdc_UUIDs,params.GDC_token,true)
    bams_ch = download.out.bam_file
    bam_size = download.out.bamsize
    aluindex(bams_ch,bam_size,params.refseq_file,params.expression_file,params.genome_file,params.snp_file,params.genome_regions_index_EI,params.regions_EI)
    aluindex.out.results.view()
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


