#!/usr/bin/env bash

# Set the environment
[ -e /opt/lofarsoft/lofarinit.sh ] && source /opt/lofarsoft/lofarinit.sh

# Run the requested command
if [ -z "$*" ]; then
  exec /bin/bash
else
  exec "$@"
fi
