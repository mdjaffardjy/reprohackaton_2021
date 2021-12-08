#!/bin/sh
docker pull mjdaff/r-base:4.0.3
docker build src/.docker_modules/r-base/4.0.3 -t 'mdjaff/r-base:4.0.3'
docker push mdjaff/r-base:4.0.3
#docker buildx build --platform linux/amd64,linux/arm64 -t "mdjaff/r-base:4.0.3" --push src/.docker_modules/r-base/4.0.3
