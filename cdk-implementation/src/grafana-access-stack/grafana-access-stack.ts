import { CfnOutput, Stack, StackProps } from 'aws-cdk-lib';
import { Effect, Policy, PolicyStatement, User } from 'aws-cdk-lib/aws-iam';
import { Construct } from 'constructs';

interface GrafanaAccessStackProps extends StackProps {
  database: string;
}

export class GrafanaAccessStack extends Stack {
  constructor(scope: Construct, id: string, props: GrafanaAccessStackProps) {
    super(scope, id, props);

    // Grafana Cloud User
    const user = new User(this, 'grafana-user', {});
    user.attachInlinePolicy(
      new Policy(this, 'read-timestream', {
        statements: [
          new PolicyStatement({
            effect: Effect.ALLOW,
            actions: [
              'timestream:DescribeDatabase',
              'timestream:ListTables',
              'timestream:Select',
              'timestream:PrepareQuery',
              'timestream:DescribeTable'
            ],
            resources: ['arn:aws:timestream:' + props.env.region + ':' + props.env.account + ':database/' + props.database'],
          }),
        ],
      })
    );

    new CfnOutput(this, 'USER_FOR_GRAFANA_CLOUD', {
      value: user.userName,
    });
  }
}
