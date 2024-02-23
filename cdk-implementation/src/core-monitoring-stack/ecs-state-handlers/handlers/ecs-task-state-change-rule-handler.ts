/* eslint-disable import/no-extraneous-dependencies */
import {
  EventBridgeClient,
  PutEventsCommand,
} from '@aws-sdk/client-eventbridge';
import {
  QueryCommand,
  QueryCommandOutput,
  TimestreamQueryClient,
} from '@aws-sdk/client-timestream-query';

import { Handler } from 'aws-lambda';
import {
  CPUAggregatedValues,
  MemAggregatedValues,
  TaskAggregatedData,
} from '../internal-models/timestream-aggregated-data-model';

// Constants
const ECS_TO_BATCH_CPU = 1024;
const MS_TO_NS = 1e6;

const timestreamQueryClient = new TimestreamQueryClient({});
const eventBridge = new EventBridgeClient({});
/**
 * The labels that we want to consider in the aggregated metrics when querying the database.
 */
const attributeLabels = JSON.parse(process.env.AGGREGATED_LABELS || '[]');

/**
 * Utility function to assert necessary variables exist,
 * and throw an appropriate error if they aren't.
 */
function validateEnvVar(varName: string): string {
  const variable = process.env[varName];
  if (!variable) {
    throw new Error(
      `${varName} is not defined. Please add it to the task environment variables.`
    );
  }
  return variable;
}

function findId(
  taskEnvs:
    | {
        name: string;
        value: string;
      }[]
    | undefined,
  taskIdLabel: string[] | string
) {
  if (taskEnvs === undefined) {
    console.log('No Task ID Label found in the event');
    return;
  }

  let taskIdValue;
  // check if taskIdLabels is a string or an array
  if (typeof taskIdLabel === 'string') {
    taskIdValue = taskEnvs.find((envVar) => envVar.name === taskIdLabel)?.value;
  } else {
    // if it is an array, then we need to find the taskIdLabel in the array
    taskIdValue = taskEnvs.find((envVar) =>
      taskIdLabel.includes(envVar.name)
    )?.value;
  }

  if (taskIdValue === undefined) {
    console.log(`No Task ID Label ${taskIdLabel} found in the event`);
    return;
  }
  return taskIdValue;
}

/**
 * The ID of the task that we want to consider in the aggregated metrics when querying the database.
 */
const taskIdLabel = validateEnvVar('TASK_ID_LABEL');

/**
 * The database from where we query the aggregated metrics.
 */
const databaseName = validateEnvVar('DATABASE_NAME');

/**
 * The event bus name where we send the aggregated metrics.
 */
const eventBusName = validateEnvVar('EVENT_BUS_NAME');

const handler: Handler<ECSTaskEventBusEvent> = async (event) => {
  console.log('Event received', JSON.stringify(event));

  const { lastStatus, startedAt, stoppedAt, overrides, cpu } = event.detail;
  const containerOverrides = overrides.containerOverrides[0];
  const taskCpus = Number.parseInt(cpu);

  // Extract required information from the event
  // const stateChange = event.detail.lastStatus;
  // const timestamp = new Date(event.time).getTime();
  // const startedAt = new Date(event.detail.startedAt).getTime();
  // const stoppedAt = new Date(event.detail.stoppedAt).getTime();
  // const durationMs = stoppedAt - startedAt;
  // const containerOverrides = event.detail.overrides.containerOverrides[0];
  // const taskCpus = Number.parseInt(event.detail.cpu);

  const taskIdValue = findId(containerOverrides.environment, taskIdLabel);

  // const taskIdValue = containerOverrides.environment?.find(
  //   (envVar) => envVar.name === taskIdLabel
  // )?.value;

  if (taskIdValue === undefined) {
    console.log(`No Task ID Label(s) ${taskIdLabel} found in the event`);
    return;
  }

  // If the event state is done, then we can query the database to get the aggregated metrics
  if (lastStatus === 'STOPPED') {
    if (stoppedAt === undefined) {
      console.log(
        'StoppedAt time is missing. Skip event. This should not happen.'
      );
      return;
    }

    const durationMs =
      new Date(stoppedAt).getTime() - new Date(startedAt).getTime();
    // take the stoppedAt time as the event timestamp. This is the time when the task stopped. If the values is missing use the current time.

    const filteredAttributeLabels =
      await filterAttributesNotAvailableInDatabase(
        attributeLabels,
        'docker_container_cpu'
      );
    const cpuAggregated = await getCpuAggregates(
      taskIdLabel,
      taskIdValue,
      durationMs,
      taskCpus,
      filteredAttributeLabels
    );
    if (cpuAggregated === undefined) {
      console.log('No CPU aggregated values found');
    }
    console.log('CPU aggregated values', JSON.stringify(cpuAggregated));

    const memAggregated = await getMemAggregates(taskIdLabel, taskIdValue);
    if (memAggregated === undefined) {
      console.log('No Mem aggregated values found');
    }
    console.log('Mem aggregated values', JSON.stringify(memAggregated));

    // Create the aggregated data object
    const aggregatedData: TaskAggregatedData = {
      jobId: taskIdValue,
      stopCode: event.detail.stopCode,
      exitCode: event.detail.containers[0].exitCode || -1,
      startedTime: new Date(startedAt).getTime(),
      stoppedTime: new Date(stoppedAt).getTime(),
      duration: durationMs,
      cpu: cpuAggregated,
      mem: memAggregated,
      labels: cpuAggregated?.labels || [],
    };

    console.log('Aggregated data', JSON.stringify(aggregatedData));
    // Send this to the event bus
    await sendToEventBus(aggregatedData);
  }
};

export {
  extractValueFromResponse,
  filterAttributesNotAvailableInDatabase,
  getCpuAggregates,
  handler
};

function extractValueFromResponse(
  response: QueryCommandOutput,
  measureName: string,
  position: number
): number {
  const scalarValueStr = response.Rows?.find(
    (value) => value.Data?.[0]?.ScalarValue === measureName
  )?.Data?.[position]?.ScalarValue;
  return scalarValueStr ? parseFloat(scalarValueStr) : -1;
}

function extractLabelValuesFromResponse(
  response: QueryCommandOutput,
  measureName: string,
  columnIndex: number,
  additionalDimensions: string[]
): { value: number | null; labels: { name: string; value: string }[] } {
  let value: number | null = null;
  const labels: { name: string; value: string }[] = [];

  if (response.Rows === undefined) {
    return { value, labels };
  }

  for (const row of response.Rows) {
    const data = row.Data;
    if (data === undefined) {
      continue;
    }
    if (data[0].ScalarValue === measureName) {
      value = parseFloat(data[columnIndex].ScalarValue as string);

      // Extract values for additional dimensions and store them in labels
      for (let i = 0; i < additionalDimensions.length; i++) {
        labels.push({
          name: additionalDimensions[i],
          value: data[i + 1].ScalarValue as string,
        });
      }
      break;
    }
  }

  return { value, labels };
}

function getPosition(position: number, additionalDimensions: string[]): number {
  return position + additionalDimensions.length;
}

async function getCpuAggregates(
  taskIdName: string,
  taskIdValue: string,
  durationMs: number,
  taskCpus: number,
  additionalDimensions: string[]
): Promise<
  | (CPUAggregatedValues & { labels: { name: string; value: string }[] })
  | undefined
> {
  const escapedDimensions = additionalDimensions.map(
    (dimension) => `"${dimension}"`
  );
  const dimensions = escapedDimensions.join(', ');
  const query = `
  SELECT
    measure_name, ${dimensions},
    COALESCE( MIN(measure_value::bigint),MIN(measure_value::double) ) AS min,
    COALESCE( AVG(measure_value::bigint),AVG(measure_value::double) ) AS avg,
    COALESCE( MAX(measure_value::bigint),MAX(measure_value::double) ) AS max
  FROM
    "${databaseName}"."docker_container_cpu"
  WHERE
    ${taskIdName} = '${taskIdValue}'
  GROUP BY
    measure_name${additionalDimensions.length > 0 ? ', ' : ''}${dimensions}`;

  const params = { QueryString: query };
  console.log('query', JSON.stringify(query));

  try {
    const command = new QueryCommand(params);
    const response = await timestreamQueryClient.send(command);
    console.log('ResponseJson:' + JSON.stringify(response));

    const USAGE_PERCENT = 'usage_percent';
    const USAGE_TOTAL = 'usage_total';
    const THROTTLING_TIME = 'throttling_throttled_time';

    const AVG = getPosition(2, additionalDimensions);
    const MAX = getPosition(3, additionalDimensions);

    const cpuAggregated: CPUAggregatedValues & {
      labels: { name: string; value: string }[];
    } = {
      assignedCount: taskCpus / ECS_TO_BATCH_CPU,
      maxPercentage: extractValueFromResponse(response, USAGE_PERCENT, MAX),
      avgPercentage: extractValueFromResponse(response, USAGE_PERCENT, AVG),
      timeReserved: (taskCpus / ECS_TO_BATCH_CPU) * durationMs * MS_TO_NS,
      timeUsed: extractValueFromResponse(response, USAGE_TOTAL, MAX),
      timeThrottled: extractValueFromResponse(response, THROTTLING_TIME, MAX),
      labels: extractLabelValuesFromResponse(
        response,
        USAGE_PERCENT,
        MAX,
        additionalDimensions
      ).labels,
    };

    return cpuAggregated;
  } catch (error) {
    console.error(error);
  } finally {
    console.log('done');
  }
  return undefined;
}

async function getMemAggregates(
  taskIdName: string,
  taskIdValue: string
): Promise<MemAggregatedValues | undefined> {
  const query = `
  SELECT
    measure_name,
    COALESCE( MIN(measure_value::bigint),MIN(measure_value::double) ) AS min,
    COALESCE( AVG(measure_value::bigint),AVG(measure_value::double) ) AS avg,
    COALESCE( MAX(measure_value::bigint),MAX(measure_value::double) ) AS max
  FROM
    "${databaseName}"."docker_container_mem"
  WHERE
    ${taskIdName} = '${taskIdValue}'
  GROUP BY
    measure_name`;

  const params = { QueryString: query };
  console.log('query', JSON.stringify(query));

  try {
    const command = new QueryCommand(params);
    const response = await timestreamQueryClient.send(command);
    console.log('ResponseJson:' + JSON.stringify(response));

    const LIMIT = 'limit';
    const USAGE = 'usage';
    const MAX_USAGE = 'max_usage';

    const MAX = 3;
    const AVG = 2;

    const max_usage = extractValueFromResponse(response, MAX_USAGE, MAX);
    const limit = extractValueFromResponse(response, LIMIT, MAX);
    const memAggregated: MemAggregatedValues = {
      limit: limit,
      avg: Math.ceil(extractValueFromResponse(response, USAGE, AVG)),
      max: max_usage,
    };
    return memAggregated;
  } catch (error) {
    console.error(error);
  } finally {
    console.log('done');
  }
  return undefined;
}

// Describe the timestream table schema and remove the attributeLabels not available in the table schema.
async function filterAttributesNotAvailableInDatabase(
  initialLabels: string[],
  tableName: string
): Promise<string[]> {
  const query = `DESCRIBE "${databaseName}"."${tableName}"`;
  const params = { QueryString: query };
  console.log('query', JSON.stringify(query));

  try {
    const command = new QueryCommand(params);
    const response = await timestreamQueryClient.send(command);
    console.log('ResponseJson:' + JSON.stringify(response));

    if (response.Rows) {
      const columnNames = response.Rows.map(
        (row) => row.Data && row.Data[0]?.ScalarValue
      );

      const filteredAttributeLabels = initialLabels.filter((value) =>
        columnNames.includes(value)
      );

      // Get elements which are not included in the filteredAttributeLabels
      const excludedElements = initialLabels.filter(
        (value) => !filteredAttributeLabels.includes(value)
      );

      // Log out the excluded elements
      console.log('Excluded Elements: ', excludedElements);

      return filteredAttributeLabels;
    }
  } catch (error) {
    console.error(error);
  } finally {
    console.log('done');
  }
  return [];
}

async function sendToEventBus(aggregatedData: TaskAggregatedData) {
  const input = {
    Entries: [
      {
        Source: 'metrics.aggregator',
        Resources: [],
        DetailType: 'TaskAggregatedData',
        Detail: JSON.stringify(aggregatedData),
        EventBusName: eventBusName,
      },
    ],
  };
  const command = new PutEventsCommand(input);
  const response = await eventBridge.send(command);
  console.log('Send to EventBus ResponseJson:' + JSON.stringify(response));
}
