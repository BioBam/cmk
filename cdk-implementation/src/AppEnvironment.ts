import { Environment } from 'aws-cdk-lib';

export interface MonitoringAgentEnv {
  /**
   * @type {number}
   * @memberof MonitoringAgent
   * @description The number of days you want to retain monitoring agent log events in the specified log group.
   */
  logRetentionDays: number;
  /**
   * @type {number}
   * @memberof MonitoringAgent
   * @default 994 (default for amazon linux 2)
   * @allowedValues 994 (default for amazon linux 2), 497 (default for amazon linux 1) , etc
   * @description
   * The root id of the ami image where the agent will finally run.
   * This is used to set the user of the container to telegraf:rootId
   * This is required because the telegraf container needs to access the docker.sock file to read statistics,
   * which is owned by root
   */
  rootId: number;
}

export interface ClusterEnv {
  vpcId: string;
  names: string[];
  sgs: string[];
}

export interface MetricsProperties {
  /**
   * Unique ENV variable for every task (AWS_BATCH_JOB_ID given by default in AWS)
   * Must be available in the task change event to identify the task
   * In multi-try with different batch task id use the original id to see all tries chart.
   */
  taskIdLabel: string;

  /**
   * Labels to keep from the running docker container labels or environment variables.
   */
  labelsToKeep: string[];

  /**
   * Exclude containers by name from the monitoring.
   * The values can have * as wildcard.
   */
  excludeContainerNames: string[];
}

export interface TimestreamProperties {
  /**
   * The name of the database where the monitoring system will store the metrics.
   */
  database: string;
  /**
   * The number of hours to retain high-resolution data in the memory store.
   */
  memoryRetentionHours: number;
  /**
   * The number of days to retain high-resolution data in the magnetic store.
   */
  magneticRetentionDays: number;
}

export interface ApiProperties {
  /**
   * Api key value to restrict access if no identity pool is used.
   */
  key: string;
}

export interface CoreMonitoringEnv {
  /**
   * @type {Environment}
   * @memberof AppEnvironment
   * @description
   * The environment where the monitoring system will be deployed
   */
  env: Environment;

  /**
   * @type {ClusterEnv}
   * @memberof AppEnvironment
   * @description
   * The cluster contains information of the ecs clusters where the monitoring agent and will run.
   */
  cluster: ClusterEnv;

  /**
   * Tags to add to all the newly created resources during the deployment of the monitoring system.
   */
  infrastructureTags: { [key: string]: string };

  /**
   * The timestream database used for the monitoring system and some table properties.
   */
  timestream: TimestreamProperties;

  /**
   * Properties for the api access.
   */
  api: ApiProperties;

  /**
   * @type {MonitoringAgentEnv}
   * @memberof AppEnvironment
   * @description
   * The monitoring agent is a container that runs on each ec2 instance and collects metrics from the host and the docker containers running on the host.
   */
  agent: MonitoringAgentEnv;

  /**
   * The attributes used for the aggregatted statistics
   */
  // attributeLabels: string[];

  metrics: MetricsProperties;
}



export interface AppEnvironment extends CoreMonitoringEnv {}
