import { CfnOutput } from 'aws-cdk-lib';
import { LambdaIntegration, RestApi } from 'aws-cdk-lib/aws-apigateway';
import { Effect, Policy, PolicyStatement } from 'aws-cdk-lib/aws-iam';
import { Architecture, Runtime } from 'aws-cdk-lib/aws-lambda';
import { NodejsFunction } from 'aws-cdk-lib/aws-lambda-nodejs';
import { RetentionDays } from 'aws-cdk-lib/aws-logs';
import { Construct } from 'constructs';

export interface StatsApiProps {
  databaseName: string;
  apiKey: string;
}

export class StatsApiConstruct extends Construct {
  constructor(scope: Construct, id: string, props: StatsApiProps) {
    super(scope, id);
    const apiKeyValue = props.apiKey;

    const getTaskExecutionStatsHandler = new NodejsFunction(
      this,
      'get-task-execution-stats',
      {
        runtime: Runtime.NODEJS_16_X,
        bundling: {
          nodeModules: ['@aws-sdk/client-timestream-query'],
        },
        architecture: Architecture.ARM_64,
        entry: `${__dirname}/handlers/get-task-execution-stats-handler.ts`,
        logRetention: RetentionDays.ONE_DAY,
        environment: { TSDB_NAME: props.databaseName },
      }
    );

    getTaskExecutionStatsHandler.role?.attachInlinePolicy(
      new Policy(this, 'read-timestream', {
        statements: [
          new PolicyStatement({
            effect: Effect.ALLOW,
            actions: ['timestream:*'],
            resources: ['*'],
          }),
        ],
      })
    );

    // Define the API
    const api = new RestApi(this, 'task-stats-api');
    api.addApiKey('api-key', {
      value: apiKeyValue, // set an api key to avoid open access.
    });

    const tasks = api.root.addResource('tasks');
    const task = tasks.addResource('{task_id}');

    const getStatsIntegration = new LambdaIntegration(
      getTaskExecutionStatsHandler
    );
    task.addMethod('GET', getStatsIntegration, {
      apiKeyRequired: true,
    });

    new CfnOutput(this, 'API_URL', {
      value: api.url,
    });
    new CfnOutput(this, 'API_KEY', {
      value: apiKeyValue,
    });
  }
}
