#!/bin/bash

# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# Adjusts container user to match host UID/GID (creates a non‑root user if
# needed) and executes the command via gosu.

set -e

# Get UID and GID from environment (set by start_container.sh or .env)
HOST_UID=${HOST_UID:-0}
HOST_GID=${HOST_GID:-$HOST_UID}   # default GID to UID if not provided

# Move libraries to third-party directory (common to both root and non‑root)
# eigen3
rm -rf /workspace/third-party/eigen3
mkdir -p /workspace/third-party/eigen3/include
cp -r /opt/include/eigen3/* /workspace/third-party/eigen3/include
# spdlog
rm -rf /workspace/third-party/spdlog
mkdir -p /workspace/third-party/spdlog/include/spdlog
cp -r /usr/include/spdlog/* /workspace/third-party/spdlog/include/spdlog
# fmt
rm -rf /workspace/third-party/fmt
mkdir -p /workspace/third-party/fmt/include/fmt
cp -r /usr/include/fmt/* /workspace/third-party/fmt/include/fmt

if [ "$HOST_UID" != "0" ]; then
    EXISTING_USER=$(getent passwd $HOST_UID | cut -d: -f1)
    if [ -z "$EXISTING_USER" ]; then
        # Create group if needed
        if ! getent group $HOST_GID >/dev/null; then
            groupadd -g $HOST_GID devgroup
            GROUP_NAME=devgroup
        else
            GROUP_NAME=$(getent group $HOST_GID | cut -d: -f1)
        fi

        # Create the user and add to sudo group
        useradd -m -u $HOST_UID -g $GROUP_NAME -G sudo -s /bin/bash devuser
        USER_NAME=devuser

        # Give passwordless sudo rights (no password prompt)
        echo "$USER_NAME ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/$USER_NAME
        chmod 440 /etc/sudoers.d/$USER_NAME
    else
        USER_NAME=$EXISTING_USER
        # If the user already exists,
        # ensure they are in the sudo group and have sudoers file
        usermod -aG sudo $USER_NAME
        if [ ! -f /etc/sudoers.d/$USER_NAME ]; then
            echo "$USER_NAME ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/$USER_NAME
            chmod 440 /etc/sudoers.d/$USER_NAME
        fi
    fi

    # Ensure home directory ownership
    if [ -d "/home/$USER_NAME" ]; then
        chown $HOST_UID:$HOST_GID "/home/$USER_NAME"
    fi

    # Run command as devuser (which can now use sudo)
    exec gosu $USER_NAME "$@"
else
    # Run as root directly
    exec "$@"
fi
