#!/bin/bash
# Set the root directory where your folders are located
root_dir="."

# Find all the subdirectories under the root directory
find "$root_dir" -mindepth 1 -type d -print0 | while IFS= read -r -d '' dir; do
    # Extract the U and V values from the directory name
    u_value=$(echo "$dir" | grep -oP 'U=-?[0-9]+(\.[0-9]+)?' | cut -d'=' -f2)
    v_value=$(echo "$dir" | grep -oP 'V=-?[0-9]+(\.[0-9]+)?' | cut -d'=' -f2)

    # Create the new directory structure if it doesn't exist
    new_dir="$root_dir/U=$u_value/V=$v_value"
    mkdir -p "$new_dir"
    echo $new_dir
    # Move the subdirectory to the new location
    mv "$dir"/* "$new_dir"/
done
