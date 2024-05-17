## Examples of Labels

A list of labels we find useful to add to the metrics of a production compute environment processing scientific workflows. All or part of the labels can be useful depending on the use case.

<!-- A table with labels, descriptions and use cases -->

| Label          | Description                                                | Use Case                                                         |
|----------------|------------------------------------------------------------|------------------------------------------------------------------|
| `task_id`      | Unique identifier for a specific task within a workflow.   | To track the performance or status of individual tasks.          |
| `workflow_id`  | Identifier for the entire workflow process.                | To aggregate metrics across all tasks in a particular workflow.  |
| `user_id`      | ID associated with the user who initiated the workflow.    | For usage tracking and auditing purposes.                        |
| `node_id`      | Identifier for the compute node executing a task.          | To identify load distribution and potential node-specific issues.|
| `project_id`   | Project-specific identifier for resource segregation.      | To facilitate billing and resource allocation per project.       |
| `environment`  | Label specifying the environment (dev, test, prod, etc.).  | To differentiate metrics from various stages of deployment, or compute environments.      |
| `data_set_id`  | ID of the dataset being processed by a task or workflow.   | To trace data usage and lineage throughout the processing chain. |
| `job_type`     | Type or classification of job (e.g., CPU-bound, I/O-heavy).| To analyze the efficiency of different types of workloads.       |
| `billing_code`      | Code used for attributing costs to a particular budget.         | To allocate costs and manage budgeting for projects or departments that use the compute resources.    |
| `partition`         | Computational partition or queue the job was submitted to.      | For analyzing system utilization and load across different computational partitions or queues.        |
| `tool_id`           | Identifier for the software or tool being executed.              | To analyze performance or usage patterns of specific tools within workflows.                          |
| `tool_version`      | Version number of the tool.                                     | For tracking outcomes and behaviors associated with different versions of a tool.                     |
| `container_image`   | The name or identifier of the container image used.             | To identify the execution context regarding the base environment of a task.                           |
| `container_tag`     | Tag associated with the container image, often a version tag.   | To specify and track which version of the container is utilized.                                      |
| `execution_engine`  | Label denoting the engine (e.g., Docker, Singularity) used.     | To compare metrics across different container runtimes or execution engines.                          |



Using these labels can help in organizing and filtering metrics meaningfully. This categorization enables quick identification of issues, assesses system performance, and assists in making informed decisions regarding scaling, maintenance, and optimization of the compute environment and resources assigned to each task.
