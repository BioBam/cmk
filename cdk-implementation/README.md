# CMK Cloud Application

## Configuration
The env folder contains an example environment configuration file.
To configure a deployment, copy the file and adjust the parameters to your environment.
The configuration file path is set in the `CONFIG` environment variable for the execution.
The value must correspond to the path to a configuration .yaml file.
e.g. CONFIG=./env/sandbox.yaml
The Yaml configuration file follows the structure of the class `AppEnvironment` defined in `src/AppEnvironments.ts`.

Example call the CDK using the sandbox configuration.
```bash
CONFIG=./env/sandbox.yaml cdk diff
```

Example call cdk using a different "production" configuration.
```bash
CONFIG=./env/prod.yaml cdk diff
```

## Projen
This project is configured as a projen project.
Projen takes care of all the dependencies and other parts of configuration.
The projen is configured in the .projenrc file.

- [Projen API Reference](https://projen.io/api/API.html)
- [Projen GitHub Repository](https://github.com/projen/projen)


# Example Workflows

The folder [bioinformatics-workflows](../bioinformatics-workflows) contains some example workflows.

Check the [evaluation readme](../bioinformatics-workflows/EVALUATION_ENVIRONMENT_README.md) for details of a use case.

