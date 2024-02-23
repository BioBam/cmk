const { awscdk } = require('projen');
const project = new awscdk.AwsCdkTypeScriptApp({
  cdkVersion: '2.72.1',
  defaultReleaseBranch: 'main',
  name: 'cmk',
  authorName: 'Robert Nica',
  authorEmail: 'rnica@biobam.com',
  packageManager: 'npm',
  runtime: awscdk.LambdaRuntime.NODEJS_18_X,
  bundlingOptions: {
    // list of node modules to exclude from the bundle
    externals: ['aws-sdk'],
    sourcemap: true,
  },

  eslint: true,
  prettier: true,
  prettierOptions: {
    settings: {
      singleQuote: true,
    },
  },

  // deps: [],                /* Runtime dependencies of this module. */
  // description: undefined,  /* The description is just a string that helps people understand the purpose of the package. */
  devDeps: [
    '@aws-sdk/client-dynamodb',
    '@aws-sdk/lib-dynamodb',
    '@aws-sdk/client-ec2',
    '@aws-sdk/client-ecs',
    '@aws-sdk/client-s3',
    '@aws-sdk/client-eventbridge',
    '@types/aws-lambda',
    '@aws-sdk/client-timestream-query',
    '@aws-sdk/client-timestream-write',
    '@types/js-yaml',
    'js-yaml',
  ],
  /* Build dependencies for this module. */
  // packageName: undefined,  /* The "name" in package.json. */
});
project.addGitIgnore('.history/');
project.addGitIgnore('.DS_Store');
project.addGitIgnore('.nextflow/');

// add eslint rules, remove the import issue in ts functions.
project.eslint.addRules({
  'import/no-extraneous-dependencies': ['error', { devDependencies: true }],
});

project.synth();
