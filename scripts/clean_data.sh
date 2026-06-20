#!/bin/bash
# Clean data directory: keep only original.txt in each test case folder
#
# Usage:
#   ./clean_data.sh          # Interactive (lists what will be deleted first)
#   ./clean_data.sh --dry-run # Show what would be deleted
#   ./clean_data.sh --force  # Delete immediately without prompting

DATA_DIR="data"
DRY_RUN=false
FORCE=false

# Parse arguments
for arg in "$@"; do
    case $arg in
        --dry-run)
            DRY_RUN=true
            ;;
        --force)
            FORCE=true
            ;;
    esac
done

# Check if data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo "Error: $DATA_DIR directory not found"
    exit 1
fi

# Counters
total_dirs=0
files_to_delete=()
errors=0

# Iterate through each subdirectory
for dir in "$DATA_DIR"/*/; do
    # Skip if not a directory
    [ -d "$dir" ] || continue
    
    total_dirs=$((total_dirs + 1))
    
    # Find all files except original.txt
    for file in "$dir"*; do
        # Skip if not a file
        [ -f "$file" ] || continue
        
        # Skip original.txt
        basename=$(basename "$file")
        if [ "$basename" = "original.txt" ]; then
            continue
        fi
        
        files_to_delete+=("$file")
    done
done

if [ ${#files_to_delete[@]} -eq 0 ]; then
    echo "No files to delete. Data directories are already clean."
    exit 0
fi

echo "Found ${#files_to_delete[@]} files to delete across $total_dirs directories."
echo ""

if [ "$DRY_RUN" = true ]; then
    echo "DRY RUN - Files that would be deleted:"
    for f in "${files_to_delete[@]}"; do
        echo "  $f"
    done
    echo ""
    echo "Run without --dry-run to actually delete these files."
    exit 0
fi

if [ "$FORCE" = false ]; then
    echo "Files to be deleted:"
    for f in "${files_to_delete[@]}"; do
        echo "  $f"
    done
    echo ""
    read -p "Delete ${#files_to_delete[@]} files? [y/N] " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Cancelled."
        exit 1
    fi
fi

# Delete files
deleted=0
for f in "${files_to_delete[@]}"; do
    rm -f "$f"
    if [ $? -eq 0 ]; then
        deleted=$((deleted + 1))
    else
        errors=$((errors + 1))
    fi
done

echo "Done!"
echo "Files deleted: $deleted"
echo "Errors: $errors"
