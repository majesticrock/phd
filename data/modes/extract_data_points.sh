#!/bin/bash

# Specify the directory where your folders are located
directory="./square/dos_900/T=0.0"

# Use the find command to locate the folders with names matching the pattern
find "$directory" -type d -name 'U=-2.0_V=*' | while read -r folder; do
    # Use awk to extract the second number
    second_number=$(echo "$folder" | awk -F'[_=]' '{print $6}')
    echo "$second_number"
done
