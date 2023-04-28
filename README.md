
# Cloud Monitoring Kit (CMK)

Cloud Monitoring Kit (CMK) is a fine-grained monitoring system optimized for executing scientific algorithms and workflows, specifically in bioinformatics, on cloud platforms such as Amazon Web Services (AWS). It offers comprehensive insights into resource utilization for individual tasks and groups of tasks, with a focus on cost-effectiveness and integration with other systems.

## Table of Contents

- [Cloud Monitoring Kit (CMK)](#cloud-monitoring-kit-cmk)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Features](#features)
  - [License](#license)

## Introduction

Compute environments are crucial for running scientific algorithms and workflows, but effective monitoring remains a challenge, particularly in cloud environments, due to insufficient metrics from cloud providers and application-specific tools. This repository contains the reference implementation of the monitoring system described in the research conducted by our team.

## Features

- Serverless, cost-effective implementation on AWS
- Monitoring tasks executed on AWS Batch, a popular backend platform for scientific workflow management systems (WMS)
- Grouping tasks by labels for aggregated statistics
- Access metrics manually or programmatically for integration with other systems

## License

This project is licensed under the Apache License Version 2.0. See the [LICENSE](LICENSE) file for details.


<br>

> **_NOTE:_**  This project is under active development. The code and additional resources will be released periodically. Keep an eye on this space for updates and enhancements.