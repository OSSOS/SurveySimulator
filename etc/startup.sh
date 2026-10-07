#!/bin/bash
# skaha launch passthrough: skaha supplies the per-session command (JupyterLab for a
# notebook session, xterm for a desktop-app session, or a script for a headless job)
# as the container CMD; we set up the environment and exec it unchanged.
set -a
. /etc/profile
exec "$@"
