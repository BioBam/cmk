/* eslint-disable import/no-extraneous-dependencies */
import {
  TimestreamWriteClient,
  WriteRecordsCommand,
} from '@aws-sdk/client-timestream-write';

import { EventBridgeEvent, Handler } from 'aws-lambda';
import { TaskAggregatedData } from '../../ecs-state-handlers/internal-models/timestream-aggregated-data-model';

const timestreamWriteClient = new TimestreamWriteClient({});

/**
 * The database from where we query the aggregated metrics.
 */
const databaseName = process.env.DATABASE_NAME;

/**
 * The table from where we query the aggregated metrics.
 */
const tableName = process.env.TABLE_NAME;

{
  if (databaseName === undefined) {
    throw new Error('DATABASE_NAME is not defined');
  }
  if (tableName === undefined) {
    throw new Error('TABLE_NAME is not defined');
  }
}

const handler: Handler<
  EventBridgeEvent<'TaskAggregatedData', TaskAggregatedData>
> = async (event) => {
  console.log('Event received', JSON.stringify(event));

  // Extract required information from the event
  const aggregatedData: TaskAggregatedData = event.detail;
  console.log('Aggregated data', aggregatedData);

  //create the dimensions using the selected labels
  const dimensions: Array<Dimension> = [];
  dimensions.push({
    Name: 'jobId',
    Value: aggregatedData.jobId,
  });
  aggregatedData.labels.forEach((label) => {
    // if (attributeLabels.includes(label.name)) {
    dimensions.push({
      Name: label.name,
      Value: label.value,
    });
    // }
  });

  // Add the exitCode and stopCode to the dimensions.
  // This is useful to filter out the metrics when querying the database.
  dimensions.push({
    Name: 'exitCode',
    Value: aggregatedData.exitCode.toString(),
  });
  dimensions.push({
    Name: 'stopCode',
    Value: aggregatedData.stopCode.toString(),
  });
  // add startTime and StoppedTime to the dimensions
  dimensions.push({
    Name: 'startedTime',
    Value: aggregatedData.startedTime.toString(),
  });
  dimensions.push({
    Name: 'stoppedTime',
    Value: aggregatedData.stoppedTime.toString(),
  });

  const measures: Measure[] = createMeasureArray(aggregatedData);

  // Insert the aggregated metrics into the database
  await insertMetrics(
    databaseName,
    tableName,
    dimensions,
    measures,
    aggregatedData.stoppedTime
  );

  // return the result of the push operation
};

export { handler };

type Dimension = {
  Name: string;
  Value: string;
};

type Measure = {
  name: string;
  value: string;
  type: string;
};

async function insertMetrics(
  dbName: string,
  tName: string,
  dimensions: Array<Dimension>,
  measures: Measure[],
  stoppedTime: number
) {
  // Set up the WriteRecordsCommand with your data
  const writeRecordsCommand = new WriteRecordsCommand({
    DatabaseName: dbName,
    TableName: tName,
    Records: measures.map((m) => ({
      Dimensions: dimensions,
      MeasureName: m.name,
      MeasureValue: m.value,
      MeasureValueType: m.type,
      Time: stoppedTime.toString(),
      TimeUnit: 'MILLISECONDS',
    })),
  });

  try {
    // Send the command to insert metrics in Timestream
    await timestreamWriteClient.send(writeRecordsCommand);
    console.log('Successfully inserted metric');
  } catch (error) {
    console.error('Error inserting metric:', error);
  }
}

function createMeasureArray(aggregatedData: TaskAggregatedData): Measure[] {
  const m: Measure[] = [];
  // duration measure
  addBigIntMeasure(m, 'durationMs', aggregatedData.duration);
  // memory measures
  const memMetrics = aggregatedData.mem;
  if (memMetrics !== undefined) {
    addBigIntMeasure(m, 'memoryMax', memMetrics.max);
    addBigIntMeasure(m, 'memoryAvg', memMetrics.avg);
    addBigIntMeasure(m, 'memoryLimit', memMetrics.limit);
  }
  const cpuMetrics = aggregatedData.cpu;
  if (cpuMetrics !== undefined) {
    addBigIntMeasure(m, 'cpuAssignedCount', cpuMetrics.assignedCount);
    addDoubleMeasure(m, 'cpuAvgPercentage', cpuMetrics.avgPercentage);
    addDoubleMeasure(m, 'cpuMaxPercentage', cpuMetrics.maxPercentage);
    addBigIntMeasure(m, 'cpuTimeReserved', cpuMetrics.timeReserved);
    addBigIntMeasure(m, 'cpuTimeThrottled', cpuMetrics.timeThrottled);
    addBigIntMeasure(m, 'cpuTimeUsed', cpuMetrics.timeUsed);
  }
  return m;
}

function addBigIntMeasure(measures: Measure[], name: string, value: number) {
  if (value != -1) {
    const intValue = Math.ceil(value);
    measures.push({
      name: name,
      value: intValue.toString(),
      type: 'BIGINT',
    });
  }
}

function addDoubleMeasure(measures: Measure[], name: string, value: number) {
  if (value != -1) {
    measures.push({
      name: name,
      value: value.toString(),
      type: 'DOUBLE',
    });
  }
}
