#!/bin/bash


tables=('docker' 'docker_container_blkio' 'docker_container_cpu' 'docker_container_mem' 'docker_container_net' 'docker_container_status' 'docker_data' 'docker_devicemapper' 'docker_metadata')

for t in "${tables[@]}"
do
    echo "Removing $t"
    aws timestream-write delete-table --database-name metricsDB --table-name $t
done