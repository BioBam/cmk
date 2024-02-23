export interface StatsValues {
  max: number;
  avg: number;
  min: number;
}

export interface TaskExecutionStats {
  cpu: {
    usagePercent: StatsValues;
    usageTotal: StatsValues;
    usageSystem: StatsValues;
    usageInUsermode: StatsValues;
    usageInKernelmode: StatsValues;
  };

  mem: {
    maxUsage: StatsValues;
    usage: StatsValues;
    limit: StatsValues;
    swap: StatsValues;
  };

  blkio: {
    // TODO: Find out what info is useful and how to use.
  };

  net: {
    // TODO: Find out what info is useful and how to use.
  };
}
