#!/bin/bash

# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# Adjusts container user to match host UID/GID (creates a non‑root user if
# needed) and executes the command via gosu.

set -e

# Get UID and GID from environment (set by start_container.sh or .env)
HOST_UID=${HOST_UID:-0}
HOST_GID=${HOST_GID:-$HOST_UID}   # default GID to UID if not provided

if [ "$HOST_UID" != "0" ]; then
    # Adjust devuser's UID to match the host user
    if [ "$(id -u devuser)" != "$HOST_UID" ]; then
        usermod -u $HOST_UID devuser
    fi
    # Adjust devgroup's GID to match the host group
    if [ "$(id -g devuser)" != "$HOST_GID" ]; then
        groupmod -g $HOST_GID devgroup
    fi
    # Fix home directory ownership after UID/GID change
    chown -R $HOST_UID:$HOST_GID /home/devuser
    exec gosu devuser "$@"
else
    exec "$@"
fi
