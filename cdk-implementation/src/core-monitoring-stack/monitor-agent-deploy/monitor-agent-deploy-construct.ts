import {
  Ec2Service,
  Ec2TaskDefinition,
  ICluster,
  PlacementConstraint,
} from 'aws-cdk-lib/aws-ecs';
import { Construct } from 'constructs';

export interface MonitoringAgentDeployProps {
  readonly cluster: ICluster;
  taskDefinition: Ec2TaskDefinition;
  serviceName: string;
}

/**
 * This construct deploys the monitoring agent as a service in the cluster,
 *  so that it gets executed on every instance running tasks in the cluster.
 *
 * @export
 * @class MonitoringAgentDeployConstruct
 * @extends {Construct}
 */
export class MonitoringAgentDeployConstruct extends Construct {
  constructor(scope: Construct, id: string, props: MonitoringAgentDeployProps) {
    super(scope, id);

    // Instantiate an Amazon ECS Service
    const telegrafDaemonService = new Ec2Service(this, 'MonitoringService', {
      cluster: props.cluster,
      taskDefinition: props.taskDefinition,
      daemon: true,
      serviceName: props.serviceName,
      placementConstraints: [PlacementConstraint.distinctInstances()],
    });
    // Instances that run someting else other than the telegraf service.
    // Note: This is a workaround for batch not scaling down to 0 issue.
    telegrafDaemonService.addPlacementConstraints(
      PlacementConstraint.memberOf('task:group != service:' + props.serviceName)
    );
  }
}
