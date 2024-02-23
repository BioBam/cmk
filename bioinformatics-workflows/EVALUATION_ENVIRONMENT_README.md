# Evaluation environment

## Table of contents

- [Evaluation environment](#evaluation-environment)
  - [Table of contents](#table-of-contents)
  - [Set up the core compute environment in AWS](#set-up-the-core-compute-environment-in-aws)
    - [Remove the configuration from AWS](#remove-the-configuration-from-aws)
  - [Set up the monitoring system in AWS](#set-up-the-monitoring-system-in-aws)
    - [Environment configuration](#environment-configuration)
      - [agent](#agent)
      - [cluster](#cluster)
      - [tags](#tags)
      - [databaseName](#databasename)
      - [attributeLabels](#attributelabels)
    - [Configure the telegraf agent](#configure-the-telegraf-agent)
      - [General Parameters](#general-parameters)
      - [Metrics, attributes, and filtering](#metrics-attributes-and-filtering)
      - [Metrics Storage](#metrics-storage)
    - [Summary of metrics](#summary-of-metrics)
    - [Deploy the CDK stacks](#deploy-the-cdk-stacks)
  - [WMS](#wms)
    - [miniwdl](#miniwdl)


## Software Requirements

- [aws cdk](https://docs.aws.amazon.com/cdk/v2/guide/getting_started.html) for deploying the CDK implementation in your AWS account
- [nextflow](https://www.nextflow.io/docs/latest/getstarted.html) for running the nextflow workflow
- [miniwdl](https://github.com/chanzuckerberg/miniwdl) for running the WDL workflow
- [miniwdl-aws](https://github.com/miniwdl-ext/miniwdl-aws) add aws batch support for miniwdl
- [vscode](https://code.visualstudio.com/) or another code editor for adjusting the configuration files


## Set up the core compute environment in AWS

To set up the AWS Core Environment we use the template "Genomics Workflow Core" from the [AWS Open Data - Genomics Workflows](https://docs.opendata.aws/genomics-workflows/quick-start.html#step-1-core-environment) Quick Start page. 
The Core Environment contains the following services needed for a running genomics workflows on the AWS Batch service:
- S3 Bucket
- EC2 Launch Templates
- IAM Roles
- AWS Batch Compute Environments
- AWS Batch Job Queues

For the configuration of the Core Environment template we use the default AWS Virtual Private Cloud (VPC) and skip the "Step 0: Amazon VPC" on the website, which creates a new VPC and would be recommended in a production environment.

Click on the "Launch Stack" button to start the CloudFormation template in our AWS account.

For the the CloudFormation parameters, we set the following:

- "VpcID" we select the default VPC from the drop-down menu.
- "NumberOfSubnets" we set to 2 and select two of the SubnetIds from the default VPC in the "SubnetIds" parameter.
- "Namespace" we set to "runworkflows" which is required later in next stacks.
- "CreateEFS" to true, to create an EFS file system for the shared file system used by miniwdl.

We leave the remaining parameters with the default values.

Following the next steps we check the acknowledge "IAM CAPABILITIES" and click on the "Create Stack" button to start the CloudFormation template.

The CloudFormation service takes a few minutes to complete the deployment.

### Remove the configuration from AWS

To remove the configuration from AWS at the end of this tutorial, we delete the stack "genomics-workflow-core" from the AWS CloudFormation service console.
This can be done after the resources created in the next steps are removed, otherwise, the CloudFormation stack deletion will fail.

## Set up the monitoring system in AWS

The monitoring system is a Cloud Development Kit (CDK) application which can be deployed on top of the "AWS Core Environment" that we deployed above.
The components of the Monitoring System can be seen in the Fig. 1. 

1. Clone the monitoring repository
2. Create a file `env/sandbox.yaml` using the template in `env/example.yaml`
3. Deploy the CDK stacks.

### Environment configuration

This is the configuration for our use case. We can edit the configuration in the `env/sandbox.yaml` file.

<details>
  <summary><code>example.yaml</code></summary>

```yaml

# Account parameters - includes the AWS region and account number
env:
  region: us-east-1  # AWS Region where resources will be deployed
  account: "012345678900"  # Your AWS account number as a string

# Clusters configuration - specifies the VPC and Security Groups to apply the monitoring system
cluster:
  vpcId: "vpc-1234567890"  # ID of your VPC
  sgs:
    - "sg-f0123459"  # IDs of the security groups attached to the clusters
  names:
    - "spot-runworkflows_Batch_1234abcd-efgh-ijkl-mn12-345abc678def"  # Name of the first ECS cluster
    # Add more ECS clusters if necessary

# Time Series Database configuration - sets up database name and retention policies
timestream:
  database: "metricsDB"  # Name of the Timestream database
  memoryRetentionHours: 1  # Duration in hours to keep the data in memory store
  magneticRetentionDays: 1000  # Duration in days to keep the data in magnetic store

# Monitoring agent configuration - defines log retention period and labels for metrics collection
agent:
  logRetentionDays: 14  # Number of days to retain logs
  rootId: 994  # Root identifier for the monitoring agent
  labels:  # Labels applied to the monitoring agent
    C_TOOL: "MonitoringAgent"

# Metrics collection configuration - specifies task ID label and list of labels to retain
metrics:
  taskIdLabel: "AWS_BATCH_JOB_ID"  # Label used to identify task ID of the batch job
  labelsToKeep:  # List of labels that should be retained when collecting metrics
    - "host"
    - "container_image"
    - "C_TOOL"
    - "C_DATASET"
    - "C_WMS"
  excludeContainerNames:  # List of container names to exclude from monitoring
    - "ecs-agent" # exclude the AWS ecs agent
    - "ecs-coremonitoringstackmonitoringagenttask*" # exclude the telegraf monitoring agent
```

</details>

### Configure the telegraf agent

The configuration for the telegraf agent is created automatically based on the parameters configured in the env yaml file.

The file `src/core-monitoring-stack/monitor-agent-task/container/telegraf.conf` is the configuration file for the agent. It is updated based on the yaml file during `cdk synth` or `cdk deploy`.

<details>
  <summary><code>telegraf.conf</code></summary>

```TOML
# ~~ Generated by the monitoring CDK app code. To modify, edit the CONFIG .yaml file and run "cdk synth"
[agent]
  interval = "10s"
  round_interval = true
  metric_batch_size = 1000
  metric_buffer_limit = 10000
  collection_jitter = "0s"
  flush_interval = "10s"
  flush_jitter = "0s"
  precision = "0s"
  hostname = "${INSTANCE_ID}"
  omit_hostname = false

[[inputs.docker]]
  endpoint = "unix:///var/run/docker.sock"
  gather_services = false
  source_tag = false
  container_name_include = []
  container_name_exclude = []
  timeout = "5s"
  perdevice = false
  perdevice_include = []
  total_include = ["cpu", "blkio", "network"]
  docker_label_include = ["AWS_BATCH_JOB_ID","C_DATASET","C_TOOL","C_WMS","host","container_image"]
  docker_label_exclude = []
  tag_env = ["AWS_BATCH_JOB_ID","C_DATASET","C_TOOL","C_WMS","host","container_image"]

[[outputs.timestream]]
  region = "${AWS_REGION}" 
  database_name = "metricsDB"
  describe_database_on_start = false
  mapping_mode = "multi-table"
  create_table_if_not_exists = true
  create_table_magnetic_store_retention_period_in_days = 7
  create_table_memory_store_retention_period_in_hours = 1
  max_write_go_routines = 25
```
</details>

To change the telegraf configuration besides the parameters in the yaml config file, edit `baseConfig` in the file [telegraf-utils.ts](../src/core-monitoring-stack/monitor-agent-task/telegraf-utils.ts) used as a base template.

#### General Parameters

Configure `[agent]` section settings for things like how often the metrics are collected. We set the instance name as the hostname for all the metrics to add the instanceID as a label to all task metrics. This is useful to then map the tasks to the instances without additional metrics plugins. 

#### Metrics, attributes, and filtering

The section `[[inputs.docker]]` contains the configuration for connecting to the docker socks and reading metrics about the containers running the tasks.
In this case we use the telegraf [docker inputs plugin](https://github.com/influxdata/telegraf/blob/master/plugins/inputs/docker/README.md) and use some parameters to specify which attributes we want to capture and which ones to filter out. We capture the `AWS_BATCH_JOB_ID` and any attributes starting with `C_` and filter out the rest. We capture these from the container tags and also from the container env. variables, as different workflows can use different ways to label the tasks.

#### Metrics Storage

We use the telegraf [outputs timestream plugin](https://github.com/influxdata/telegraf/tree/master/plugins/outputs/timestream) to output the metrics into Amazon Timestream. The `[[outputs.timestream]]` section configures this plugin.

### Aggregation of metrics

The metrics captured by the agent are summarized when the task ends. This is done in the aggregation [lambda function handler](../src/core-monitoring-stack/ecs-state-handlers/handlers/ecs-task-state-change-rule-handler.ts).

### Deploy the CDK stacks

```bash
# We can deploy all stacks using the --all parameter
CONFIG=./env/sandbox.yaml cdk deploy --all

# Or we can deploy only the Core and Grafana User by passing these stacks as parameters
CONFIG=./env/sandbox.yaml cdk deploy core-monitoring-stack monitoring-service
```


## WMS

### miniwdl

After [installing miniwdl and miniwdl-aws](./miniwdl/README.md), we can run the workflow with the miniwdl wms on the aws batch backend. We use the env variable (experimental until release) to set an environmental variable with the attribute that we need for the monitoring.
wdl has currently no support for resource labeling like nextflow does so in this case we use the environment variable to add the custom attributes to the metrics.

The workflow file is defined in the `miniwdl` folder. The input parameters are passed as command line arguments and the files have been uploaded and are available in S3.

```batch
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

### Nextflow

[Nextflow](https://www.nextflow.io/) is a DSL workflow manager.

#### Training material and documentation

- [Nextflow tutorials](https://nf-co.re/usage/nextflow)
- [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html)

#### Community-developed workflows in Nextflow

The [nf-core projects](https://nf-co.re/) hosts community curated Nextflow pipelines.

### Labeling the tasks for the monitoring system

The monitoring system needs to have the custom labels in the containier environment variables.

In nextflow workflows we can use the `containerOptions` parameter to set these. This can be added in each process, or at the top of the workflow, or in the general configuration file, depending on how the selected labels change.

In the `.nf` workflow we have added the labels in the process.

```
containerOptions "--env C_TOOL=XXXX --env C_DATASET=${params.ref} --env C_WMS=nextflow"
```

#### Running the proof of concept Nextflow pipeline

1. Install `nextflow` by following the instructions [here](https://www.nextflow.io/)
2. Optionally move the `nextflow` directory to a directory in your path
   i.e.  `echo $PATH` and then `mv nextflow usr/bin`
3. Navigate to the directory where you saved `example.nf`
4. Run the pipeline with the following command:

      ```
      nextflow run nextflow/example.nf \
            -with-docker quay.io/nextflow/rnaseq-nf:v1.0 \
            --left $PWD/test_data/reads_1.fq.gz \
            --right $PWD/test_data/reads_2.fq.gz \
            --ref $PWD/test_data/transcriptome.fa
      ```

5. Enable the use Docker containers adding the following option to the above command line: `-with-docker quay.io/nextflow/rnaseq-nf:v1.0`

Note: input files requires the use of absolute paths

#### Workflow with channels

The workflow has been modified to run for a list of datasets. Instead of a pair of reads it can be executed on a list of paired reads.

```
nextflow run example_w_channels.nf \
      -c nextflow-aws.config \
      -with-docker docker.io/robertbio/wf \
      --outdir $WF_DATA_BUCKET/channel/results_nextflow/ \
      -bucket-dir $WF_DATA_BUCKET/channel/nextflow_bucket_dir/
``````
