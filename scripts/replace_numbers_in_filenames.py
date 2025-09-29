import os
import sys
import re

def rename_files_in_folder(folder_path):
    """
    Renames files in a given folder by replacing numbers in their names
    according to a predefined mapping.
    """
    # Define the mapping from old numbers to new numbers
    replacements = {
        '72': '28', '278': '100', '434': '160', '521': '184',
        '691': '252', '936': '322', '949': '328', '1011': '350',
        '1034': '365', '1223': '426', '1351': '474', '1405': '490',
        '1514': '526', '1829': '621', '2017': '690', '2119': '722',
        '2192': '747', '2312': '790', '2454': '839', '2676': '924',
        '2700': '935', '2769': '954', '2799': '963', '2891': '993',
        '3048': '1048', '3222': '1118', '3233': '1122', '3295': '1147',
        '3338': '1160'
    }

    # 1. Validate that the folder path exists
    if not os.path.isdir(folder_path):
        print(f"❌ Error: Folder not found at '{folder_path}'")
        return

    print(f"🔍 Scanning files in '{folder_path}'...")
    renamed_count = 0

    # 2. Iterate over each file in the directory
    for original_filename in os.listdir(folder_path):
        new_filename = original_filename

        # 3. Apply all replacement rules to the current filename
        for old_num, new_num in replacements.items():
            pattern = r'\b' + old_num + r'\b'
            new_filename = re.sub(pattern, new_num, new_filename)

        # 4. If the filename was changed, rename the file
        if new_filename != original_filename:
            old_path = os.path.join(folder_path, original_filename)
            new_path = os.path.join(folder_path, new_filename)

            # Safety check to avoid overwriting another file
            if os.path.exists(new_path):
                print(f"⚠️  Skipping rename: '{new_filename}' already exists.")
                continue

            # Perform the rename
            os.rename(old_path, new_path)
            print(f"✅ Renamed: '{original_filename}'  ->  '{new_filename}'")
            renamed_count += 1

    print(f"\n✨ Done. Renamed {renamed_count} file(s).")


if __name__ == "__main__":
    # Ensure a folder path is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python rename_files.py <folder_path>")
        sys.exit(1)

    folder_to_process = sys.argv[1]
    rename_files_in_folder(folder_to_process)