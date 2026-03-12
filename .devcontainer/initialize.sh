#!/bin/bash
# This script runs before the container starts.
# It writes the host user's UID and GID into a .env file,
# which docker-compose will automatically pick up.

ENV_FILE="$(dirname "$0")/../docker/.env"

# Write (or overwrite) the .env file with the current user's IDs
echo "HOST_UID=$(id -u)" > "$ENV_FILE"
echo "HOST_GID=$(id -g)" >> "$ENV_FILE"
