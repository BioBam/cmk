import {
  QueryCommand,
  QueryCommandInput,
  TimestreamQueryClient,
} from '@aws-sdk/client-timestream-query';
import { StringUtils } from '../../../packages/utils/string-utils';
// import { Ec2InstanceChangeEventConstruct } from './ec2-instance-change-event-construct';

// const ec2Utils = new Ec2Utils();
const tsqc = new TimestreamQueryClient({});

let tsDBName = process.env.TSDB_NAME!;

export const handler = async (event: any, context: any, callback: any) => {
  console.log('Handle stats api call. env TSDB: ' + tsDBName);
  StringUtils.printEvent(event);
  const taskId = event.pathParameters.task_id; //task_id

  const params: QueryCommandInput = {
    // QueryString: 'SELECT * FROM ' + tsDBName + 'WHERE AWS_BATCH_JOB_ID = ' + taskId,
    QueryString:
      'SELECT * FROM "' +
      tsDBName +
      '"."docker_container_cpu" WHERE AWS_BATCH_JOB_ID = \'' +
      taskId +
      "' ORDER BY time DESC LIMIT 1",
  };
  const command = new QueryCommand(params);

  // async/await.
  try {
    const data = await tsqc.send(command);
    console.log('Received query result');
    StringUtils.printEvent(data);
    return data;
    // process data.
  } catch (error) {
    console.log('Received query result');
    StringUtils.printEvent(error);
    callback(Error('Failed to query.'));
  } finally {
    // finally.
  }

  // // Save information in the database
  // if ( event.detail.state === 'running') {
  //   // 1. GetDetails of instance
  //   const instance = await ec2Utils.getInstanceDetails(event.detail['instance-id']);
  //   // 2. Save in database
  //   await store.saveInstanceDetails(instance);
  // }
};

// options for return:
// https://stackoverflow.com/questions/54626183/whats-the-right-way-to-return-from-an-aws-lambda-function-in-node-js
