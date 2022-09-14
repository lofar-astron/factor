#!/bin/bash
set -e

COMMUNITY="astronrd"
REPOSITORY="factor"
BRANCH=$(git branch --show-current)
TAG=$(git describe --tags)

cd "$(dirname "${0}")"

docker build -t "${COMMUNITY}/${REPOSITORY}:${TAG}" .
docker tag "${COMMUNITY}/${REPOSITORY}:${TAG}" "${COMMUNITY}/${REPOSITORY}:${BRANCH}"

docker push "${COMMUNITY}/${REPOSITORY}:${TAG}"
docker push "${COMMUNITY}/${REPOSITORY}:${BRANCH}"
