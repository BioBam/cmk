/**
 * This interface represents the aggregated values for a container's CPU usage.
 * It includes some CPU-specific metrics.
 */
export interface CPUAggregatedValues {
  /**
   * The number of cpus assigned to this container.
   */
  assignedCount: number;

  /**
   * The maximum percentage of CPU used. This can be more than 100% as 1 CPU is 100%.
   */
  maxPercentage: number;

  /**
   * The average percentage of CPU used. This can be more than 100% as 1 CPU is 100%.
   */
  avgPercentage: number;

  /**
   * This refers to the total amount of a specific resource used by a container.
   * It can include CPU usage in nanoseconds,
   * and block I/O usage in bytes. The usage_total metric is calculated since the container was started.
   */
  timeUsed: number;

  /**
   * This refers to the total cpu time that this container could have used.
   */
  timeReserved: number;

  /**
   * In the context of Docker statistics,
   * throttling_throttled_time refers to the total amount of time (in nanoseconds)
   * that a container has been throttled due to exceeding its CPU usage limit.
   */
  timeThrottled: number;

  /**
   * Dimensions of the task.
   */
  labels: { name: string; value: string }[];
}

/**
 * This interface represents the aggregated values for a container's memory usage.
 */
export interface MemAggregatedValues {
  /**
   * The memory limit set to this container.
   */
  limit: number;
  /**
   * The average amount of memory in bytes used by a container.
   */
  avg: number;
  /**
   * The maximum amount of memory in bytes used by a container.
   */
  max: number;
}

/**
 * This interface represents the aggregated values for a container's CPU and memory usage.
 * This is the details content of the event sent to EventBridge.
 */
export interface TaskAggregatedData {
  /**
   * The Id of the job in AWS Batch.
   */
  jobId: string;

  /**
   * The stop code message of the task in AWS Batch.
   */
  stopCode: string;

  /**
   * The exit code of the container.
   */
  exitCode: number;

  /**
   * The labels of the container.
   */
  labels: { name: string; value: string }[];

  /**
   * This refers to the total amount of time (in milliseconds) that a container has been running for.
   */
  duration: number;

  /**
   * The timestamp when the container started.
   */

  startedTime: number;

  /**
   * The timestamp when the container stopped.
   */
  stoppedTime: number;

  /**
   * The aggregated values for the CPU usage of a container.
   */
  cpu: CPUAggregatedValues | undefined;

  /**
   * The aggregated values for the memory usage of a container.
   */
  mem: MemAggregatedValues | undefined;
}

export interface TaskAggregatedDataEvent {
  Source: string;
  Resources: string[];
  DetailType: string;
  Detail: string;
  EventBusName: string;
}

export class TaskAggregatedDataEventBuilder {
  static source: string = 'metrics.aggregator';
  static detailTypeString: string = 'TaskAggregatedData';
  static build(
    detailsData: TaskAggregatedData,
    eventBusName: string
  ): TaskAggregatedDataEvent {
    return {
      Source: this.source,
      Resources: [],
      DetailType: this.detailTypeString,
      Detail: JSON.stringify(detailsData),
      EventBusName: eventBusName,
    };
  }
}
