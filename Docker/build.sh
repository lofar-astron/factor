#!/bin/bash
set -e

COMMUNITY="astronrd"
REPOSITORY="factor"
TAG=$(git describe --tags)

cd $(dirname ${0})

docker build -t ${COMMUNITY}/${REPOSITORY}:${TAG} .
docker tag ${COMMUNITY}/${REPOSITORY}:${TAG} ${COMMUNITY}/${REPOSITORY}:latest

docker push ${COMMUNITY}/${REPOSITORY}:${TAG}
docker push ${COMMUNITY}/${REPOSITORY}:latest
