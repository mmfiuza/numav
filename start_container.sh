# Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#!/bin/bash

# Add all private keys from a directory to ssh-agent
add_ssh_keys() {
    local key_dir="${1:-$HOME/.ssh}"  # Default to ~/.ssh if no argument
    
    # Check if directory exists
    if [ ! -d "$key_dir" ]; then
        echo "Error: Directory $key_dir does not exist" >&2
        return 1
    fi
    
    # Check if ssh-agent is running
    if [ -z "$SSH_AUTH_SOCK" ]; then
        echo "Error: ssh-agent is not running" >&2
        return 1
    fi
    
    # Loop through all files in directory
    for key in "$key_dir"/*; do
        # Skip if not a regular file
        [ -f "$key" ] || continue
        
        # Skip files that are definitely not private keys
        if [[ "$key" == *.pub ]] || \
           [[ "$key" == *known_hosts* ]] || \
           [[ "$key" == *authorized_keys* ]] || \
           [[ "$key" == *config* ]]; then
            continue
        fi
        
        # Skip if file is world-readable (private keys should be protected)
        if [ "$(stat -c %a "$key" | cut -c 3)" -ge 4 ]; then
            echo "Warning: Skipping $key - file is world-readable" >&2
            continue
        fi
        
        # Try to add the key
        if ssh-add "$key" 2>/dev/null; then
            echo "Added: $(basename "$key")"
        else
            echo "Skipped: $(basename "$key") (not a valid private key or requires passphrase)"
        fi
    done
}

add_ssh_keys "$1" # Usage: add_ssh_keys /path/to/keys

docker compose up -d

docker compose exec numav bash
