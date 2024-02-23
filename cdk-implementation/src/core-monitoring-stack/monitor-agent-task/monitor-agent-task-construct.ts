import * as path from 'path';
import {
  ContainerImage,
  Ec2TaskDefinition,
  LogDriver,
  NetworkMode,
} from 'aws-cdk-lib/aws-ecs';
import { Effect, PolicyStatement } from 'aws-cdk-lib/aws-iam';
import { Construct } from 'constructs';
import { TelegrafUtils } from './telegraf-utils';
import { MonitoringAgentEnv, TimestreamProperties } from '../../AppEnvironment';

export interface MonitorAgentTaskProperties {
  excludeContainerNames: string[];
  agentEnv: MonitoringAgentEnv;
  taskLabels: string[];
  timestream: TimestreamProperties;
}

/**
 * This construct creates the monitoring agent docker and pushes it into ECR.
 * After this it creates a task definition in ECS that can be used to run the agent.
 *
 * @export
 * @class TaskMonitoringAgentEcsTaskDefinitionConstruct
 * @extends {Construct}
 */

export class MonitorAgentTaskConstruct extends Construct {
  public taskDefinition: Ec2TaskDefinition;
  constructor(scope: Construct, id: string, props: MonitorAgentTaskProperties) {
    super(scope, id);

    const { agentEnv } = props;

    if (agentEnv.logRetentionDays < 1) {
      throw new Error('logRetentionDays should be a positive number');
    }

    this.taskDefinition = new Ec2TaskDefinition(this, 'MonitoringAgent', {
      networkMode: NetworkMode.HOST,
    });

    this.taskDefinition.addToTaskRolePolicy(
      new PolicyStatement({
        effect: Effect.ALLOW,
        actions: ['timestream:*'],
        resources: ['*'],
      })
    );

    // Create the telegraf.conf based on cmk config
    TelegrafUtils.createConfigFile(
      path.join(__dirname, 'container', 'telegraf.conf'),
      {
        timestream: props.timestream,
        labels: props.taskLabels,
        excludeContainerNames: props.excludeContainerNames,
      }
    );
    
    // Create the telegraf container image
    const containerImage = ContainerImage.fromAsset(
      path.join(__dirname, 'container'),
      {}
    );

    const containerDefinition = this.taskDefinition.addContainer(
      'DefaultContainer',
      {
        image: containerImage,
        memoryLimitMiB: 128,
        user: 'telegraf:' + agentEnv.rootId,
        command: [
          '/bin/bash',
          '-c',
          'export INSTANCE_ID=$(curl http://169.254.169.254/latest/meta-data/instance-id); export AWS_REGION=$(curl http://169.254.169.254/latest/meta-data/placement/region); telegraf',
        ],
        logging: LogDriver.awsLogs({
          streamPrefix: 'MonitoringAgent',
          logRetention: agentEnv.logRetentionDays,
        }),
        environment: {},
      }
    );

    // mount docker.sock as a volume inside the task
    this.taskDefinition.addVolume({
      name: 'DockerSocksAccess',
      host: {
        sourcePath: '/var/run/docker.sock',
      },
    });
    containerDefinition.addMountPoints({
      containerPath: '/var/run/docker.sock',
      sourceVolume: 'DockerSocksAccess',
      readOnly: true,
    });
  }
}
