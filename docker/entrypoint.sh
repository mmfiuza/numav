#!/bin/bash
set -e

# Get UID and GID from environment (set by start_container.sh or .env)
HOST_UID=${HOST_UID:-0}
HOST_GID=${HOST_GID:-$HOST_UID}   # default GID to UID if not provided

# If we need to switch to a non-root user
if [ "$HOST_UID" != "0" ]; then
    # Check if a user with this UID already exists
    EXISTING_USER=$(getent passwd $HOST_UID | cut -d: -f1)
    if [ -z "$EXISTING_USER" ]; then
        # Create a group with the given GID if it doesn't exist
        if ! getent group $HOST_GID >/dev/null; then
            groupadd -g $HOST_GID devgroup
            GROUP_NAME=devgroup
        else
            GROUP_NAME=$(getent group $HOST_GID | cut -d: -f1)
        fi

        # Create the user
        useradd -m -u $HOST_UID -g $GROUP_NAME -s /bin/bash devuser
        USER_NAME=devuser
    else
        USER_NAME=$EXISTING_USER
    fi

    # Ensure the home directory is owned by the user
    if [ -d "/home/$USER_NAME" ]; then
        chown $HOST_UID:$HOST_GID "/home/$USER_NAME"
    fi

    # Drop privileges and execute the command (e.g., bash)
    exec gosu $USER_NAME "$@"
else
    # Run as root if UID=0
    exec "$@"
fi