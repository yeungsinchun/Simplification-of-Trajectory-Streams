import kagglehub
import os
import shutil

# Download latest version
path = kagglehub.dataset_download("arashnic/tdriver")

print("Path to dataset files:", path)

# Check what's in the downloaded directory
print("\nContents of download:")
for item in os.listdir(path):
    full_path = os.path.join(path, item)
    if os.path.isdir(full_path):
        print(f"  [DIR] {item}/")
        # Show first few files
        subfiles = os.listdir(full_path)[:5]
        for sf in subfiles:
            print(f"         {sf}")
        if len(os.listdir(full_path)) > 5:
            print(f"         ... and {len(os.listdir(full_path)) - 5} more files")
    else:
        print(f"  {item}")
