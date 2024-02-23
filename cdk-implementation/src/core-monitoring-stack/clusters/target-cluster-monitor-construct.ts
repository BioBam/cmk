import { Aws } from 'aws-cdk-lib';
import { SecurityGroup, Vpc } from 'aws-cdk-lib/aws-ec2';
import { Cluster } from 'aws-cdk-lib/aws-ecs';
import { Construct } from 'constructs';
import { ClusterEnv } from '../../AppEnvironment';

export interface ClustersProps {
  clusterEnv: ClusterEnv;
}

/**
 * Utilities construct that retrieves the referenced clusters where the monitoring agent needs to be deployed.
 */
export class ClustersConstruct extends Construct {
  public clusters;
  constructor(scope: Construct, id: string, props: ClustersProps) {
    super(scope, id);

    const { clusterEnv } = props;

    const vpc = Vpc.fromLookup(this, 'vpc', {
      vpcId: clusterEnv.vpcId,
    });

    // Recover the security groups form the cluster
    const securityGroups = clusterEnv.sgs.map((sgId, index) =>
      SecurityGroup.fromSecurityGroupId(this, `SG${index}`, sgId)
    );

    this.clusters = clusterEnv.names.map((name, index) => {
      if (name === '') return undefined;
      const clusterName = name;
      // Get cluster reference by arn
      return Cluster.fromClusterAttributes(this, `cluster${index}`, {
        clusterArn: `arn:aws:ecs:${Aws.REGION}:${Aws.ACCOUNT_ID}:cluster/${clusterName}`,
        clusterName,
        vpc,
        securityGroups,
      });
    });
  }
}
