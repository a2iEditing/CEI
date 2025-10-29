#!/usr/bin/env bash

RESOURCES_LIST=${1:-"$(dirname "$(readlink -f "$0")")/resources_list.txt"}
RESOURCES_DIR=${2:-"Resources"}
MAIN_DIR=${3:-"main"}

missing=()

while IFS= read -r f; do
    file="$MAIN_DIR/$RESOURCES_DIR/$f"
    if [[ ! -e "$file" ]]; then
        missing+=("$file")
    fi
done < $RESOURCES_LIST

if (( ${#missing[@]} > 0 )); then
    echo "The following files are missing:"
    printf '  %s\n' "${missing[@]}"
    exit 1
else
    echo "INFO: All files exist."
    exit 0
fi
