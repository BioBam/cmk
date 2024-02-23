# Telegraf agent details

The telegraf agent runs on the ec2 instances where we need to read metrics of the docker containers. The agent runs inside a Docker container itself and has a configuration file `telegraf.conf`.

```
docker run -it --rm \
    --name telegraf \
    --user telegraf:$(stat -c '%g' /var/run/docker.sock) \
    -v /home/ec2-user/telegraf.conf:/etc/telegraf/telegraf.conf:ro \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -h $(curl http://169.254.169.254/latest/meta-data/instance-id) \
    telegraf
```

## Parameters
- **container name:** set to "telegraf" so we can exclude this container from the measurements in the `telegraf.config`
- **container user:** set to telegraf and group set to the same group as the `docker.sock`
- **configuration file:** is mounted from the user home, and contains the docker input plugin and the timestream output plugin configurations
- **docker server access:** the /var/run/docker.sock is mounted inside the container, which together with the container user belonging to the same group allows the telegraf docker input plugin to communicate with the instance docker service.
- **host name:** the telegraf container instance is set the hostname the same as the ec2 instance ID. this is to have the metrics add the "host" dimension which can then be searched in the DB.
- **image:** the image where telegraf is running is the official `telegraf` image. This is to simplify installation and configuration.

## Telegraf Configuration File
The file `telegraf.conf` contains all the parameters of the telegraf agent and configures how and what measurements are taken and where to leave them.
