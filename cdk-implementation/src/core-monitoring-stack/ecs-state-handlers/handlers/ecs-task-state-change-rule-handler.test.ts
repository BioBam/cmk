import { QueryCommandOutput } from '@aws-sdk/client-timestream-query';
import { extractValueFromResponse } from './ecs-task-state-change-rule-handler';

describe('extractValueFromResponse', () => {
  const response: QueryCommandOutput = {
    Rows: [
      {
        Data: [
          { ScalarValue: 'measureA' },
          { ScalarValue: '10' },
          { ScalarValue: '20' },
        ],
      },
      {
        Data: [
          { ScalarValue: 'measureB' },
          { ScalarValue: '30' },
          { ScalarValue: '40' },
        ],
      },
    ],
    QueryId: undefined,
    ColumnInfo: undefined,
    $metadata: {},
  };

  it('should extract the correct value from the response', () => {
    expect(extractValueFromResponse(response, 'measureA', 1)).toEqual(10);
    expect(extractValueFromResponse(response, 'measureB', 2)).toEqual(40);
  });

  it('should return -1 if the measure is not found in the response', () => {
    expect(extractValueFromResponse(response, 'measureC', 0)).toEqual(-1);
  });

  it('should return -1 if the position is out of bounds', () => {
    expect(extractValueFromResponse(response, 'measureA', 5)).toEqual(-1);
  });

  it('should return -1 if the response is invalid', () => {
    expect(
      extractValueFromResponse(
        {
          QueryId: undefined,
          Rows: undefined,
          ColumnInfo: undefined,
          $metadata: {},
        },
        'measureA',
        1
      )
    ).toEqual(-1);
  });
});

describe('extractValueFromResponse', () => {
  const response = {
    Rows: [
      {
        Data: [
          {
            ScalarValue: 'throttling_throttled_time',
          },
          {
            ScalarValue: '0.0',
          },
          {
            ScalarValue: '0.0',
          },
          {
            ScalarValue: '0.0',
          },
        ],
      },
      {
        Data: [
          {
            ScalarValue: 'usage_percent',
          },
          {
            ScalarValue: '0.030634',
          },
          {
            ScalarValue: '10.01889366391422',
          },
          {
            ScalarValue: '160.29128462311556',
          },
        ],
      },
    ],
    QueryId: undefined,
    ColumnInfo: undefined,
    $metadata: {},
  };
  it('should extract the correct value from the response', () => {
    expect(
      extractValueFromResponse(response, 'throttling_throttled_time', 3)
    ).toEqual(0);
    expect(extractValueFromResponse(response, 'usage_percent', 1)).toEqual(
      0.030634
    );
  });
});
