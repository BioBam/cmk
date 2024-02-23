import { Duration, Stack } from 'aws-cdk-lib';
import { EventBus, Rule } from 'aws-cdk-lib/aws-events';
import { LambdaFunction } from 'aws-cdk-lib/aws-events-targets';
import { Effect, Policy, PolicyStatement } from 'aws-cdk-lib/aws-iam';
import { Architecture, Runtime } from 'aws-cdk-lib/aws-lambda';
import { NodejsFunction } from 'aws-cdk-lib/aws-lambda-nodejs';
import { RetentionDays } from 'aws-cdk-lib/aws-logs';
import { Construct } from 'constructs';
import { MetricsProperties } from '../../AppEnvironment';

export interface EcsStateHandlersProps {
  clusterNames: string[];
  metricsProps: MetricsProperties;
  stack: Stack;
  eventBus: EventBus;
  databaseName: string;
}

export class EcsStateHandlersConstruct extends Construct {
  constructor(scope: Construct, id: string, props: EcsStateHandlersProps) {
    super(scope, id);

    // EventBridge rule to match the tasks state changes in the selected clusters where monitoring is active
    const ecsTaskChangedEventRule = new Rule(this, 'anyStateChangeRule', {
      description: 'Changes in the aws ecs task state',
      enabled: true,
      eventPattern: {
        source: ['aws.ecs'],
        detailType: ['ECS Task State Change'],
        detail: {
          clusterArn: props.clusterNames.map(
            (cluster) =>
              'arn:aws:ecs:' +
              props.stack.region +
              ':' +
              props.stack.account +
              ':cluster/' +
              cluster
          ),
        },
      },
    });

    const ecsTaskDoneHandler = new NodejsFunction(
      this,
      'ecs-task-state-change-handler',
      {
        runtime: Runtime.NODEJS_18_X,
        bundling: {
          nodeModules: [
            '@aws-sdk/client-timestream-query',
            '@aws-sdk/client-eventbridge',
          ],
        },
        architecture: Architecture.ARM_64,
        entry: `${__dirname}/handlers/ecs-task-state-change-rule-handler.ts`,
        logRetention: RetentionDays.ONE_DAY,
        timeout: Duration.seconds(300),
        memorySize: 256,
        environment: {
          DATABASE_NAME: props.databaseName,
          TASK_ID_LABEL: props.metricsProps.taskIdLabel,
          AGGREGATED_LABELS: JSON.stringify(props.metricsProps.labelsToKeep),
          EVENT_BUS_NAME: props.eventBus.eventBusName,
        },
      }
    );

    // Allow the lambda to write to the event bus
    props.eventBus.grantPutEventsTo(ecsTaskDoneHandler);

    // Allow the lambda to read from timestream
    ecsTaskDoneHandler.role?.attachInlinePolicy(
      new Policy(this, 'write-timestream', {
        statements: [
          new PolicyStatement({
            effect: Effect.ALLOW,
            actions: ['timestream:*'],
            resources: ['*'],
          }),
        ],
      })
    );

    ecsTaskChangedEventRule.addTarget(
      new LambdaFunction(ecsTaskDoneHandler, {
        // deadLetterQueue: queue, // Optional: add a dead letter queue
        maxEventAge: Duration.hours(2), // Optional: set the maxEventAge retry policy
        retryAttempts: 0, // Optional: set the max number of retry attempts
      })
    );
  }
}
