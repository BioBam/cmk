import { EventBus, Rule } from 'aws-cdk-lib/aws-events';
import { CloudWatchLogGroup } from 'aws-cdk-lib/aws-events-targets';
import { LogGroup, RetentionDays } from 'aws-cdk-lib/aws-logs';
import { Construct } from 'constructs';

export interface EventBusLogProps {
  eventBus: EventBus;
}

export class EventBusLogConstruct extends Construct {
  constructor(scope: Construct, id: string, props: EventBusLogProps) {
    super(scope, id);

    // EventBridge rule to match the metrics aggregator events and log to CloudWatch Logs
    const anyEventRule = new Rule(this, 'any-event-rule', {
      enabled: true,
      eventPattern: {
        source: ['metrics.aggregator'],
      },
      eventBus: props.eventBus,
    });

    const logGroup = new LogGroup(this, 'bus', {
      retention: RetentionDays.ONE_DAY,
    });
    anyEventRule.addTarget(new CloudWatchLogGroup(logGroup));
  }
}
