#!/bin/sh
docker pull lbmc/htseq:0.13.5
# docker build src/.docker_modules/htseq/0.13.5 -t 'lbmc/htseq:0.13.5'
# docker push lbmc/htseq:0.13.5
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/htseq:0.13.5" --push src/.docker_modules/htseq/0.13.5
