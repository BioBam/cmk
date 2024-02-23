import { Stack, StackProps } from 'aws-cdk-lib';
import { EventBus } from 'aws-cdk-lib/aws-events';
import { Construct } from 'constructs';
import { CoreMonitoringEnv } from '../AppEnvironment';
import { ClustersConstruct } from './clusters/target-cluster-monitor-construct';
import { EcsStateHandlersConstruct } from './ecs-state-handlers/ecs-state-handlers-construct';
import { EventBusAggregatesHandlersConstruct } from './event-bus-aggregates-handlers/event-bus-aggregates-handlers-construct';
import { EventBusLogConstruct } from './event-bus-aggregates-log/event-bus-log-aggregates-construct';
import { MonitoringAgentDeployConstruct } from './monitor-agent-deploy/monitor-agent-deploy-construct';
import { MonitorAgentTaskConstruct } from './monitor-agent-task/monitor-agent-task-construct';
import { StatsApiConstruct } from './stats-api/stats-api-construct';

interface MonitoringSystemProps extends StackProps {
  coreEnv: CoreMonitoringEnv;
}

/**
 * This stack defines the monitoring system that watches the tasks and has information of the execution statistics of these.
 */
export class MonitoringSystemStack extends Stack {
  constructor(scope: Construct, id: string, props: MonitoringSystemProps) {
    super(scope, id, props);
    const { coreEnv } = props;

    // Create the monitoring agent task definition
    const ecsDefinitionConstruct = new MonitorAgentTaskConstruct(
      this,
      'monitoring-agent-task',
      {
        excludeContainerNames: coreEnv.metrics.excludeContainerNames,
        agentEnv: coreEnv.agent,
        taskLabels: [
          ...coreEnv.metrics.labelsToKeep,
          coreEnv.metrics.taskIdLabel,
        ],
        timestream: coreEnv.timestream,
      }
    );

    // Retrieve the clusters where the monitoring agent needs to be deployed
    const { clusters } = new ClustersConstruct(this, 'clusters', {
      clusterEnv: coreEnv.cluster,
    });

    // Deploy the monitoring agent in each cluster
    clusters.forEach((clusterReference, index) => {
      if (clusterReference === undefined) return;
      new MonitoringAgentDeployConstruct(this, `telegrafC${index}`, {
        cluster: clusterReference,
        taskDefinition: ecsDefinitionConstruct.taskDefinition,
        serviceName: `telegrafC${index}`,
      });
    });

    // Create the event bus to send the events to
    const eventBus = new EventBus(this, 'monitoring-event-bus', {
      eventBusName: 'monitoring',
    });

    // Create the handlers for the task state change events
    new EcsStateHandlersConstruct(this, 'ecs-notification', {
      clusterNames: clusters.map(
        (cluster) => cluster?.clusterName ?? 'NO_MATCH'
      ),
      metricsProps: coreEnv.metrics,
      stack: this,
      eventBus: eventBus,
      databaseName: coreEnv.timestream.database,
    });

    // Log aggregated metrics
    new EventBusLogConstruct(this, 'event-bus-log', {
      eventBus,
    });

    // Push aggregated metrics into aggregated_metrics Timestream DB Table
    new EventBusAggregatesHandlersConstruct(this, 'event-bus-aggregates', {
      eventBus,
      timestream: coreEnv.timestream,
    });

    // Create the REST API to access metrics
    new StatsApiConstruct(this, 'api', {
      databaseName: coreEnv.timestream.database,
      apiKey: coreEnv.api.key
    });
  }
}
