import { Duration } from 'aws-cdk-lib';
import { EventBus, Rule } from 'aws-cdk-lib/aws-events';
import { LambdaFunction } from 'aws-cdk-lib/aws-events-targets';
import { Effect, Policy, PolicyStatement } from 'aws-cdk-lib/aws-iam';
import { Architecture, Runtime } from 'aws-cdk-lib/aws-lambda';
import { NodejsFunction } from 'aws-cdk-lib/aws-lambda-nodejs';
import { RetentionDays } from 'aws-cdk-lib/aws-logs';
import { CfnTable } from 'aws-cdk-lib/aws-timestream';
import { Construct } from 'constructs';
import { TimestreamProperties } from '../../AppEnvironment';
import { TaskAggregatedDataEventBuilder } from '../ecs-state-handlers/internal-models/timestream-aggregated-data-model';

export interface EventBusAggregatesHandlersProps {
  eventBus: EventBus;
  timestream: TimestreamProperties;
}

export class EventBusAggregatesHandlersConstruct extends Construct {
  constructor(
    scope: Construct,
    id: string,
    props: EventBusAggregatesHandlersProps
  ) {
    super(scope, id);

    const table = new CfnTable(this, 'aggregated_metrics', {
      tableName: 'aggregated_metrics',
      databaseName: props.timestream.database,
      retentionProperties: {
        memoryStoreRetentionPeriodInHours:
          props.timestream.memoryRetentionHours,
        magneticStoreRetentionPeriodInDays:
          props.timestream.magneticRetentionDays,
      },
    });

    // EventBridge rule to match the arrival of a new metrics aggregator event
    const anyEventRule = new Rule(this, 'any-metrics-event', {
      enabled: true,
      eventPattern: {
        source: [TaskAggregatedDataEventBuilder.source],
        detailType: [TaskAggregatedDataEventBuilder.detailTypeString],
      },
      eventBus: props.eventBus,
    });

    const pushAggregatesIntoDatabase = new NodejsFunction(
      this,
      'task-aggregates-handler',
      {
        runtime: Runtime.NODEJS_18_X,
        bundling: {
          nodeModules: ['@aws-sdk/client-timestream-write'],
        },
        architecture: Architecture.ARM_64,
        entry: `${__dirname}/handlers/push-into-database-handler.ts`,
        logRetention: RetentionDays.ONE_DAY,
        timeout: Duration.seconds(300),
        memorySize: 256,
        environment: {
          TABLE_NAME: table.tableName!,
          DATABASE_NAME: table.databaseName,
          EVENT_BUS_NAME: props.eventBus.eventBusName,
        },
      }
    );

    pushAggregatesIntoDatabase.role?.attachInlinePolicy(
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

    anyEventRule.addTarget(
      new LambdaFunction(pushAggregatesIntoDatabase, {
        maxEventAge: Duration.hours(2),
        retryAttempts: 0,
      })
    );
  }
}
