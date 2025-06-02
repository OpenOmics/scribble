#!/bin/bash

# Script to find the common parent directory of multiple file paths
# Usage: ./common_parent.sh path1 path2 path3 ...
# Returns: The longest common parent directory path

set -euo pipefail

# Function to find the common parent directory of multiple file paths
find_common_parent() {
    local paths=("$@")
    
    # Check if any paths provided
    if [[ ${#paths[@]} -eq 0 ]]; then
        echo "Error: No paths provided" >&2
        return 1
    fi
    
    # If only one path, return its directory
    if [[ ${#paths[@]} -eq 1 ]]; then
        dirname "${paths[0]}"
        return 0
    fi
    
    # Convert all paths to absolute paths and get their directories
    local abs_dirs=()
    local i
    for i in "${!paths[@]}"; do
        local path="${paths[i]}"
        
        # Get absolute path
        if [[ -e "$path" ]]; then
            # File/directory exists - use realpath for accurate resolution
            if command -v realpath >/dev/null 2>&1; then
                if [[ -f "$path" ]]; then
                    abs_dirs+=("$(dirname "$(realpath "$path")")")
                else
                    abs_dirs+=("$(realpath "$path")")
                fi
            else
                # Fallback if realpath not available
                if [[ -f "$path" ]]; then
                    abs_dirs+=("$(cd "$(dirname "$path")" && pwd)")
                else
                    abs_dirs+=("$(cd "$path" && pwd)")
                fi
            fi
        else
            # File doesn't exist - work with the path as-is
            if [[ "$path" = /* ]]; then
                # Already absolute
                if [[ "$path" =~ /[^/]+\.[^/]*$ ]]; then
                    # Looks like a file (has extension)
                    abs_dirs+=("$(dirname "$path")")
                else
                    # Assume it's a directory
                    abs_dirs+=("$path")
                fi
            else
                # Make relative path absolute
                local abs_path="$(pwd)/$path"
                if [[ "$path" =~ /[^/]+\.[^/]*$ ]]; then
                    # Looks like a file
                    abs_dirs+=("$(dirname "$abs_path")")
                else
                    # Assume it's a directory
                    abs_dirs+=("$abs_path")
                fi
            fi
        fi
    done
    
    # Split the first directory into components
    local first_dir="${abs_dirs[0]}"
    IFS='/' read -ra first_parts <<< "$first_dir"
    
    # Find common prefix
    local common_parts=()
    local part_index=0
    
    # Iterate through each part of the first directory
    for part in "${first_parts[@]}"; do
        local is_common=true
        
        # Check if this part exists in all other directories at the same position
        for dir in "${abs_dirs[@]:1}"; do
            IFS='/' read -ra dir_parts <<< "$dir"
            
            # Check if this directory has enough parts and the part matches
            if [[ ${#dir_parts[@]} -le $part_index ]] || [[ "${dir_parts[$part_index]}" != "$part" ]]; then
                is_common=false
                break
            fi
        done
        
        if [[ "$is_common" == true ]]; then
            common_parts+=("$part")
        else
            break
        fi
        
        ((part_index++))
    done
    
    # Reconstruct the common path
    local common_path=""
    for part in "${common_parts[@]}"; do
        if [[ -n "$part" ]]; then
            common_path+="/$part"
        fi
    done
    
    # Handle edge cases
    if [[ -z "$common_path" ]]; then
        common_path="/"
    fi
    
    echo "$common_path"
}

# Print usage information
usage() {
    cat << EOF
Usage: $0 <file_path1> [file_path2] [file_path3] ...

Find the common parent directory of multiple file paths.

Arguments:
    file_path1, file_path2, ...    File or directory paths to analyze

Examples:
    $0 /home/user/data/file1.txt /home/user/data/file2.txt
    $0 ./project/src/main.c ./project/include/header.h ./project/docs/readme.md
    $0 /var/log/app.log /var/cache/app.cache /var/lib/app.db

Output:
    Prints the longest common parent directory path to stdout.

Exit codes:
    0 - Success
    1 - Error (no arguments provided or other error)
EOF
}

# Main script logic
main() {
    # Check for help flag
    if [[ $# -eq 1 && ("$1" == "-h" || "$1" == "--help") ]]; then
        usage
        exit 0
    fi
    
    # Check if any arguments provided
    if [[ $# -eq 0 ]]; then
        echo "Error: No file paths provided" >&2
        echo >&2
        usage >&2
        exit 1
    fi
    
    # Validate that all arguments are provided
    for arg in "$@"; do
        if [[ -z "$arg" ]]; then
            echo "Error: Empty path argument provided" >&2
            exit 1
        fi
    done
    
    # Find and print the common parent directory
    echo $(find_common_parent "$@")
}

# Run main function with all command line arguments
main "$@"