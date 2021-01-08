#!/bin/bash
# remove exited containers:
docker ps --filter status=dead --filter status=exited -aq | xargs -r docker rm -v
# remove unused containers:
yes | docker container prune
# remove unused images:
yes | docker image prune
# remove unused volumes:
yes | docker volume prune
echo ""
df -h
