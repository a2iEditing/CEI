# CEI
A Cytoplasmic Index for Quantifying Immune-Related A-to-I RNA Editing

This repository contains both the **cloud-computing platform** and **analysis code** associated with our academic paper:

> **"TODO"**  
> Authors
> [DOI: 10.1234/example.doi.2025.001](https://doi.org/10.1234/example.doi.2025.001)


## 📑 Table of Contents
- [Overview](#-overview)
- [Repository Structure](#-repository-structure)
- [Computational Cloud-Native Pipeline](#computational-cloud---native-pipeline)
  - [Getting Started](#-getting-started---initialization)
    - [Step 1 - Initialize the Machine](#step-1---initialize-the-machine)
    - [Step 2 - Initialize the Bucket](#step-2---initialize-the-bucket)
    - [Step 3 - Initialize Resources](#step-3---initialize-resources)
  - [Run the Analysis](#run-the-analysis)
  - [Profiles](#profiles)
  - [Parameters Details](#parameters-details)
    - [AWS](#aws)
    - [GCP](#gcp)
  - [Docker Images](#docker-images)
    - [Workflow](#Workflow)
    - [Initialization](#Initialization)
- [Computation of Inverted and Tandem Repeats](#computation-of-inverted-and-tandem-repeats)

---

## 🧬 Overview

This project offers:
- A cloud-native computing platform built with Nextflow, Docker, and GCP/AWS for high-throughput RNA editing analysis
- Resources for quantifying both the **cytoplasmic editing index (CEI)** and **global Alu editing index (AEI)**
- Reproducible workflows for the analysis presented in the paper
- Example datasets and configurations for running your own analysis

---

## 📁 Repository Structure
CEI provides:
- Cloud-Native Pipeline
  - A scalable cloud-native workflow for computing:
    - The **Cytoplasmic Editing Index (CEI)**  
    - The **Global Alu Editing Index ([AEI](https://www.nature.com/articles/s41592-019-0610-9))**
    - Editing within known editing sites for [human](https://doi.org/10.1038/s41467-022-28841-4) and [mouse](https://doi.org/10.1186/gb-2014-15-1-r5)
    - Gene expression using [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html)
  - Infrastructure support for:
    - **AWS and GCP** - note that the AWS pipeline is more advanced and we recommend using this version.
    - Parallelized processing of **SRA** (FASTQ) and **TCGA** (BAM) samples.
  - Full analysis and plotting scripts for:
    - Benchmarking CEI vs global index
    - Immune-related editing response
    - Orientation-aware Alu element analysis
- Analysis
  - Python script for generating inverted and tandem repeats, including Conda `.yml` file with enviroment requirements.
  - All scripts for reproducing article figures and annotations

---
## Computational Cloud-Native Pipeline         

### 🚀 Getting Started - Initialization

Before running the workflows, you must follow the initialization sequence below.

#### Step 1 - Initialize the Machine
The machine should have at least 64GB RAM, number of cores matching your chosen `NUM_THREADS` parameter in the initialization script command (16 cores is recommended).
Storage requirements for initialization: 350GB.

- AWS - ensure that your AWS machine role is configured with the necessary policies as outlined in the [Nextflow AWS documentation](https://www.nextflow.io/docs/edge/aws.html).


Make sure your machine has updated versions of [Nextflow](https://www.nextflow.io/docs/latest/install.html) and [Docker](https://docs.docker.com/engine/install/). Make sure your user account has permission to access the Docker daemon.    
In addition, the following should be available: `wget`, `curl`, `gzip`, `awk`. Notice `zip` and  `unzip` are also required for Nextflow installation.     

Per-platform requirements:    
- AWS - download and initialize [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).
This project utilizes the AWS Batch Fargate service as the execution environment. It is essential to configure an instance role for the Batch service with the following permissions:
    - AmazonS3FullAccess
    - AmazonEC2ContainerServiceforEC2Role
    - AmazonEC2ContainerRegistryFullAccess
    - AmazonEC2FullAccess
    - AmazonECS_FullAccess
    - AWSBatchServiceRole
  - Nextflow on AWS Fargate leverages [Seqera Wave](https://seqera.io/wave). Ensure you have a Seqera (Tower) token as outlined in the [Wave documentation](https://www.nextflow.io/docs/edge/wave.html#wave-page).
  - For additional guidance, refer to the [Nextflow Fargate documentation](https://www.nextflow.io/docs/edge/aws.html#aws-fargate) and the [AWS Fargate documentation](https://docs.aws.amazon.com/batch/latest/userguide/fargate.html). 
- GCP - download and initialize [gcloud CLI](https://cloud.google.com/sdk/docs/install) (including `gsutil`), initialize [GCP credentials](https://nextflow.io/docs/latest/google.html#cloud-batch) and make sure [GCP Batch API](https://cloud.google.com/batch/docs/get-started) is enabled. Configure Google credentials and ensure the Google Cloud Batch API is enabled as described in the [Nextflow documentation](https://www.nextflow.io/docs/latest/google.html#google-cloud).

#### Step 2 - Initialize the Bucket
Create a bucket for resources and output.      
- AWS - create a bucket that your machine role and batch instance role can access according to the policies you set from [Nextflow documentation](https://www.nextflow.io/docs/edge/aws.html)
- GCP - uniform access to bucket and the following IAM permissions:

| **Principal**            | **Role**                             |
|--------------------------|--------------------------------------|
| **Editors of the project**| Storage Legacy Bucket Owner, Storage Legacy Object Owner         |
| **Owners of the project** | Storage Legacy Bucket Owner, Storage Legacy Object Owner         |
| **Service Account**       | Compute Engine Service Agent |
| **Viewers of the project**| Storage Legacy Bucket Reader, Storage Legacy Object Reader        |


#### Step 3 - Initialize Resources
Run the resource initialization script to initialize all resources for human (hg38) and mm10 (mouse):
```bash
nohup sh CloudPipeline/Init/init_main.sh <PLATFORM> <BUCKET_NAME> <NUM_THREADS>  > init.out 2> init.err &
```
| Parameter  | Position | Description | Allowed Values |
|------------|----------|-------------|----------------|
| `PLATFORM` | 1 | Which platform is used | `GCP` or `AWS` |
| `BUCKET_NAME` | 2 | Resources bucket name | Pre-created bucket name (not including `s3://` or `gs://`) |
| `NUM_THREADS` | 3 | Number of threads for STAR and salmon index generation | Positive whole numbers |           

Default for number of threads for generation of STAR and Salmon indices is 10. STAR requires at least 64G RAM for this process.      
For other organisms, you must supply all resources downloaded in pipeline and them generate STAR and Salmon indices.

### Run the Analysis
After initialization, launch the analysis workflow.     
Example for AWS - for GCP, use the files within the GCP directory.      
[Profiles](#profiles) and [parameters](#parameters-details) detailed below.        

**1. Update user parameters configuration file:**    
Change the user parameters within ``CloudPipeline/AWS/SRA_pipeline/rna_editing.user_params.config``

**2. Run:**
```bash
nohup ~/nextflow -c CloudPipeline/AWS/SRA_pipeline/rna_editing.config -bg run CloudPipeline/AWS/SRA_pipeline/rna_editing.nf -profile <SE,stranded,RL75,hg38> --run_title <RUN_TITLE> -bucket-dir <NEXTFLOW_BUCKET_WORKDIR> --srrACC_list <SRR_LIST> > log.out 2> log.err &
```
For restricted access data (dbGaP)
```bash
nohup ~/nextflow -c CloudPipeline/AWS/SRA_pipeline/rna_editing.config -bg run CloudPipeline/AWS/SRA_pipeline/rna_editing.nf -profile <SE,stranded,RL75,hg38> --run_title <RUN_TITLE> -bucket-dir <NEXTFLOW_BUCKET_WORKDIR> --srrACC_list <SRR_LIST> --NGC_file <NGC_FILE> > log.out 2> log.err &
```
Note the `-bucket-dir <NEXTFLOW_BUCKET_WORKDIR>` is only required for AWS.        

---
    
### Profiles
The following profiles are supported. Any combination of profile options can be used from these categories, but at least one profile must be used from each category. See [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles) to learn more about profiles.       

| Category  | Options | Description | Steps Affected |
|-----------|---------|-------------|----------------|
| Genome | `hg38`, `mm10` | Genome of organism - human or mouse | Alignment, expression & editing quantification |
| Library type | `SE`, `PE` | Sequencing library type - single-end or paired-end | Preprocessing, alignment, expression & editing quantification |
| Sequencing directionality | `stranded`, `unstranded` | Is the sequencing directional? | Editing quantification |
| Sequencing length | `RL50`, `RL75`, `RL100`, `RL125`, `RL150` | Read length | Preprocessing & alignment |        

Note that stranded data can also be run as unstranded, without utilizing the strand information.     
Read length should match the general read length or below, as reads shorten than wanted length by 3bp or more will be filtered out (see protocol at (TODO add article reference)[].

    
---
    
### Parameters Details
The following parameters must be provided, either directly via flags or within a configuration file:
    
#### AWS        
For using AWS to run on SRA (see [Nextflow for Fargate documentation](https://www.nextflow.io/docs/latest/aws.html#aws-fargate)):
| Parameter  | Description | Type   | Configuration File |
|------------|-------------|--------|--------------------|
| `--ecr_region` | AWS region | AWS parameter | User parameters config |
| `--process_queue` | AWS Batch queue for Fargate | AWS parameter | User parameters config |
| `--jobRole` | AWS Batch job role | AWS parameter | User parameters config |
| `--executionRole` | AWS Batch execution role | AWS parameter | User parameters config |
| `--tower_access_token` | Nextflow [Seqera access token](https://www.nextflow.io/docs/latest/wave.html) | Nextflow parameter | User parameters config |
| `--workspace_id` | Nextflow [Seqera workspace ID](https://www.nextflow.io/docs/latest/wave.html) | Nextflow parameter | User parameters config |
| `--bucket_name` | Resources bucket | AWS parameter | User parameters config |
| `--run_title` | Title of current run | Per-run parameter | General config |
| `--srrACC_list` | File with SRA accessions | Per-run parameter | General config |
| `--NGC_file` | NGC file for [restricted dbGaP access](https://www.ncbi.nlm.nih.gov/sra/docs/sra-dbGAP-cloud-download/), not required | Per-run parameter | General config |

For using AWS to run on TCGA (see [Nextflow for Fargate documentation](https://www.nextflow.io/docs/latest/aws.html#aws-fargate)):
| Parameter  | Description | Type   | Configuration File |
|------------|-------------|--------|--------------------|
| `--ecr_region` | AWS region | AWS parameter | User parameters config |
| `--process_spot_queue` | AWS Batch queue for Fargate | AWS parameter | User parameters config |
| `--jobRole` | AWS Batch job role | AWS parameter | User parameters config |
| `--executionRole` | AWS Batch execution role | AWS parameter | User parameters config |
| `--tower_access_token` | Nextflow [Seqera access token](https://www.nextflow.io/docs/latest/wave.html) | Nextflow parameter | User parameters config |
| `--workspace_id` | Nextflow [Seqera workspace ID](https://www.nextflow.io/docs/latest/wave.html) | Nextflow parameter | User parameters config |
| `--resources_bucket_name` | Resources bucket | AWS parameter | User parameters config |
| `--results_bucket_name` | Output bucket | AWS parameter | User parameters config |
| `--GDC_token` | GCD [authorization token](https://docs.gdc.cancer.gov/Data/Data_Security/Data_Security/) | User parameter | User parameters config |
| `--run_title` | Title of current run in general config | Per-run parameter | General config |
| `--gdc_UUID_list` | File with UUID accessions | Per-run parameter | General config |
     
#### GCP      
For using GCP to run on SRA:     
| Parameter  | Description | Type   | Configuration File |
|------------|-------------|--------|--------------------|
| `--project_name` | GCP project workspace name | GCP parameter | User parameters config |
| `--region` | GCP region | GCP parameter | User parameters config |
| `--bucket_name` | Resources bucket | GCP parameter | User parameters config |
| `--run_title` | Title of current run in general config | Per-run parameter | General config |
| `--srrACC_list` | File with SRA accessions | Per-run parameter | General config |

    
---
    

### Docker Images
#### Workflow
| Step  | Image | Source | 
|-------|-------|--------|
| Download SRA FASTQ files | levanonlab/sratoolkit:3.2.1 | [staphb/sratoolkit:3.2.1](https://hub.docker.com/r/staphb/sratoolkit/) |
| Preprocessing | levanonlab/fastp:0.23.4--hadf994f_2 | [quay.io/biocontainers/fastp:0.23.4--hadf994f_2](https://quay.io/repository/biocontainers/fastp?tab=tags&tag=0.23.4--hadf994f_2) |
| Expression quantification | levanonlab/salmon:1.10.2--hecfa306_0 | [quay.io/biocontainers/salmon:1.10.2--hecfa306_0](https://quay.io/repository/biocontainers/salmon?tab=tags&tag=1.10.2--hecfa306_0) |
| Alignment | levanonlab/star:2.7.10b--h9ee0642_0 | [quay.io/biocontainers/star:2.7.10b--h9ee0642_0](https://quay.io/repository/biocontainers/star?tab=tags&tag=2.7.10b--h9ee0642_0) |
| RNA editing index | levanonlab/rna-editing-index-lite:1.0 |   |
| RNA editing index & timing | levanonlab/rna-editing-index-lite:1.0.time |   |
| Per-site editing quantification | levanonlab/cmpileup:1.0 |   |
| Download TCGA BAM files | levanonlab/gdc-client:2.3.0 |  |

#### Initialization       
| Step  | Image | Source | 
|-------|-------|--------|
| Resource processing | levanonlab/bedtools:2.31.1 | [quay.io/biocontainers/bedtools:2.31.1--h13024bc_3](https://quay.io/repository/biocontainers/bedtools?tab=tags&tag=latest) |
| Resource processing | levanonlab/bigbedtobed:482--h0b57e2e_0 | [quay.io/biocontainers/ucsc-bigbedtobed:482--h0b57e2e_0](https://quay.io/repository/biocontainers/bedtools?tab=tags&tag=latest) |
| Region index generation | levanonlab/rna-editing-index-lite:1.0 | |
| Salmon index generation | levanonlab/salmon:1.10.2--hecfa306_0 | [quay.io/biocontainers/salmon:1.10.2--hecfa306_0](https://quay.io/repository/biocontainers/salmon?tab=tags&tag=1.10.2--hecfa306_0) |
| STAR index generation | levanonlab/star:2.7.10b--h9ee0642_0 | [quay.io/biocontainers/star:2.7.10b--h9ee0642_0](https://quay.io/repository/biocontainers/star?tab=tags&tag=2.7.10b--h9ee0642_0) |       


## Computation of Inverted and Tandem Repeats      
### Prepare
#### Enviroment
This script uses python 3.8 and [pybedtools](https://daler.github.io/pybedtools/index.html). You can install these or use the supplied `.yml` file.
#### Resources
Required resources are the repeat file and regions file, both BED6 format.        
Note that for repeat file the score column (4) has no use and can contain any information. The family column (5) can contain the repeat family names (for example, B1 or B2 for mouse) a single uniform value if no per-family computation is required (for example, "Alu" for human).       
For the regions file, the name (4) and score (5) columns are not of use.      
Both files can be given as plain text or gzipped.      
### Run
Example command:     
```bash
findOppositeOrientationRepeatsInRegions.py -r repeats_file.bed.gz -i regions_file.bed -o output_directory
```
The script gives many options for different output modes, filtering, restrictions and computations and all can be detailed using:     

```bash
findOppositeOrientationRepeatsInRegions.py --help
```
