#!/bin/bash
set -e

COMMUNITY="astronrd"
REPOSITORY="factor"
TAG=$(git describe --tags)

cd $(dirname ${0})

echo docker build -t ${COMMUNITY}/${REPOSITORY}:${TAG}
echo docker push ${COMMUNITY}/${REPOSITORY}:${TAG}
