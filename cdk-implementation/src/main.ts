import { readFileSync } from 'fs';

import { App, Stack, StackProps } from 'aws-cdk-lib';
import { Construct } from 'constructs';
import { load } from 'js-yaml';
import { AppEnvironment } from './AppEnvironment';
import { MonitoringSystemStack } from './core-monitoring-stack/monitoring-system-stack';
import { GrafanaAccessStack } from './grafana-access-stack/grafana-access-stack';

export class MyStack extends Stack {
  constructor(scope: Construct, id: string, props: StackProps = {}) {
    super(scope, id, props);
  }
}

const configFile = process.env.CONFIG || './env/example.yaml';
const appEnv = load(readFileSync(configFile, 'utf8')) as AppEnvironment;

const app = new App();

// The core monitoring stack
new MonitoringSystemStack(app, 'core-monitoring-stack', {
  coreEnv: appEnv,
  env: {
    account: appEnv.env.account,
    region: appEnv.env.region,
  },
});

// Creates a user to access the metrics DB from external Grafana Cloud service
new GrafanaAccessStack(app, 'grafana-access', {
  env: {
    account: appEnv.env.account,
    region: appEnv.env.region,
  },
  database: appEnv.timestream.database
});

app.synth();
