#!/bin/bash

# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# Runs on the host before the container starts, writes the host user's UID and
# GID into docker/.env so docker-compose can pass them to the container.

ENV_FILE="$(dirname "$0")/.env"

# Write (or overwrite) the .env file with the current user's IDs
echo "HOST_UID=$(id -u)" > "$ENV_FILE"
echo "HOST_GID=$(id -g)" >> "$ENV_FILE"

# Ensure ssh-agent fallback path exists
mkdir -p /tmp/no-ssh-agent
