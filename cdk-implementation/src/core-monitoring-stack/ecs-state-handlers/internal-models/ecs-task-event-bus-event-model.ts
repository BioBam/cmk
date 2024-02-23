interface Attribute {
  name: string;
  value: string;
}

interface Container {
  containerArn: string;
  exitCode?: number;
  lastStatus: string;
  name: string;
  image: string;
  imageDigest: string;
  runtimeId: string;
  taskArn: string;
  networkInterfaces: any[];
  cpu: string;
  memory: string;
}

interface ContainerOverride {
  command?: string[];
  environment?: Array<{ name: string; value: string }>;
  cpu: number;
  memory: number;
  name: string;
}

interface ECSTaskEventBusEvent {
  version: string;
  id: string;
  detailType: string;
  source: string;
  account: string;
  time: string;
  region: string;
  resources: string[];
  detail: {
    attachments: any[];
    attributes: Attribute[];
    availabilityZone: string;
    clusterArn: string;
    connectivity: string;
    connectivityAt: string;
    containerInstanceArn: string;
    containers: Container[];
    cpu: string;
    createdAt: string;
    desiredStatus: string;
    enableExecuteCommand: boolean;
    executionStoppedAt: string;
    group: string;
    launchType: string;
    lastStatus: string;
    memory: string;
    overrides: {
      containerOverrides: ContainerOverride[];
    };
    pullStartedAt: string;
    pullStoppedAt: string;
    startedAt: string;
    stoppingAt: string;
    stoppedAt: string;
    stoppedReason: string;
    stopCode: string;
    taskArn: string;
    taskDefinitionArn: string;
    updatedAt: string;
    version: number;
  };
}
