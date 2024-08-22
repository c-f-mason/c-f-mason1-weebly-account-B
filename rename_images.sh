#!/bin/bash

# Find and rename all files with 'jpg' in their name but with weird extensions
find . -type f -name '*jpg*' | while read file; do
    newfile=$(echo "$file" | sed -e 's/\?.*//')
    if [ "$file" != "$newfile" ]; then
        mv "$file" "$newfile"
        echo "Renamed: $file -> $newfile"
    fi
done

