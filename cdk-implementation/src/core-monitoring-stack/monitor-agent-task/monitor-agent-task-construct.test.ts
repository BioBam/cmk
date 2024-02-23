import { Stack } from 'aws-cdk-lib';
import { MonitorAgentTaskConstruct } from './monitor-agent-task-construct';

describe('MonitorAgentTaskDefConstruct', () => {
  test('constructor should throw if logRetentionDays is less than 1', () => {
    expect(() => {
      new MonitorAgentTaskConstruct(new Stack(), 'id', {
        agentEnv: { rootId: 994, logRetentionDays: 0 },
      });
    }).toThrow();
  });

  test('constructor should not throw if inputs are valid', () => {
    expect(() => {
      new MonitorAgentTaskConstruct(new Stack(), 'id', {
        agentEnv: { rootId: 994, logRetentionDays: 7 },
      });
    }).not.toThrow();
  });
});
