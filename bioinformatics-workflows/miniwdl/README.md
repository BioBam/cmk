# Set up clusters and workflows

## Set up AGC CLI

Follow the info on <https://aws.github.io/amazon-genomics-cli/docs/getting-started/installation/>
to set up the AGC CLI.

## Set up the cloud environment

Follow: <https://aws.github.io/amazon-genomics-cli/docs/getting-started/setup/>

(Setting up the environment via AGC CLI did not work for me so I did it using the cloudformation templates for the genomics core instead. With some work, these could be translated into CDK code).

## set up miniwdl

### install miniwdl

<https://miniwdl.readthedocs.io/en/latest/getting_started.html>

```
pip3 install miniwdl
```

### install miniwdl-aws

The AWS Batch support for the miniwdl wms. <https://github.com/miniwdl-ext/miniwdl-aws>

The experimental features like the env input support are already built in the main releases (pull request[#617](https://github.com/chanzuckerberg/miniwdl/pull/617)) but are only active in the not-yet-official WDL workflow versions 1.1 or 1.2.
The env will eventually be merged in the openwdl specification as discussed in the openwdl pull request [#504](https://github.com/openwdl/wdl/pull/504).

```
pip3 install miniwdl-aws
```

## Run using miniwdl

Environment Variables required:
```
AWS_ACCOUNT_ID = "000000000000" # The AWS account ID
AWS_REGION = us-east-1 # Region where the AWS environment is running
S3_BUCKET = s3://example-bucket # The bucket for the workflow data
```

Create an alias
```
alias awsminiwdl='miniwdl-aws-submit --workflow-queue arn:aws:batch:$AWS_REGION:$AWS_ACCOUNT_ID:job-queue/default-runworkflows --task-queue arn:aws:batch:$AWS_REGION:$AWS_ACCOUNT_ID:job-queue/default-runworkflows --no-efs --mount /mnt/efs'
awsminiwdl --s3upload $S3_BUCKET/miniwdloutput/ -w -f  workflow.wdl reads1=$S3_BUCKET/input/reads_1.fq.gz reads2=$S3_BUCKET/input/reads_2.fq.gz ref_txome=$S3_BUCKET/input/transcriptome.fa

```

Using environment variables

```
MINIWDL__AWS__WORKFLOW_QUEUE=default-runworkflows \
MINIWDL__AWS__TASK_QUEUE=default-runworkflows \
MINIWDL__AWS__FS=0 \
MINIWDL__AWS__WORKFLOW_ROLE=arn:aws:iam::$AWS_ACCOUNT_ID:role/miniwdl-workflow-role \
MINIWDL__AWS__S3UPLOAD=$S3_BUCKET/miniwdloutput/ \
miniwdl-aws-submit --mount /mnt/efs --s3upload $S3_BUCKET/miniwdloutput/ -w -f workflow.wdl \
  reads1=$S3_BUCKET/input/reads_1.fq.gz \
  reads2=$S3_BUCKET/input/reads_2.fq.gz \
  ref_txome=$S3_BUCKET/input/transcriptome.fa
```

### Details

I ran the workflow with the miniwdl wms on the aws batch backend. I used the env variable (experimental until release) to set an environmental variable with the attribute that I needed for the monitoring.

The telegraf monitoring agent needs then to listen on both, the environment variable and on tags.


##### Example run miniwdl

MINIWDL__AWS__WORKFLOW_QUEUE=default-runworkflows \
MINIWDL__AWS__TASK_QUEUE=default-runworkflows \
MINIWDL__AWS__FS=0 \
MINIWDL__AWS__WORKFLOW_ROLE=arn:aws:iam::$AWS_ACCOUNT_ID:role/miniwdl-workflow-role \
MINIWDL__AWS__S3UPLOAD=$S3_BUCKET/3step/results_miniwdl/ \
miniwdl-aws-submit --mount /mnt/efs --s3upload $S3_BUCKET/miniwdloutput/ -w -f workflow3step.wdl \
  reads1=$S3_BUCKET/simple/data_mouse/SRR937564_1.fastq.gz \
  reads2=$S3_BUCKET/simple/data_mouse/SRR937564_2.fastq.gz \
  ref=$S3_BUCKET/simple/data_mouse/mouse_ref.fa \
  output_s3_path=$S3_BUCKET/3step/results_miniwdl/

### Reference data I could use for the evaluation

#### Example1

NCBI BioProject: PRJNA419302
SRA Experiments: SRR6312174, SRR6312175, SRR6312181, SRR6312182, SRR6312187, and SRR6312190.
NCBI Nucleotide: TSA: Monilinia laxa, transcriptome shotgun assembly.

Files are available in S3 buckets directly. e.g. from the Data Access tab: s3://sra-pub-zq-6/SRR6312174/SRR6312174.sralite.1 (if downloaded from us-east-1 no charges will be incurred).

#### Example2

NCBI BioProject: PRJNA325641.  
SRA Experiments: [SRR3666079-SRR3666090](<https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP076546&o=acc_s%3Aa>)  
Reference Genome: Neisseria gonorrhoeae FA 1090 (GCA_000006845). 
