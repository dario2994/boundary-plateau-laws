import sys
import re

def replace_numbers_in_file(filename):
    """
    Replaces specific numbers in a file according to a predefined mapping.
    The replacement happens only for whole numbers.
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

    try:
        # Read the original content of the file
        with open(filename, 'r') as file:
            content = file.read()

        # Keep track of original content to see if changes were made
        original_content = content

        # Iterate through the replacement map and apply changes
        for old_num, new_num in replacements.items():
            # Use regex with word boundaries (\b) to match whole numbers only.
            # This ensures '72' is matched, but '172' or '72a' are not.
            pattern = r'\b' + old_num + r'\b'
            content = re.sub(pattern, new_num, content)

        # Write the modified content back to the file only if it has changed
        if content != original_content:
            with open(filename, 'w') as file:
                file.write(content)
            print(f"✅ Successfully processed and updated '{filename}'.")
        else:
            print(f"✅ No matching numbers found to replace in '{filename}'.")

    except FileNotFoundError:
        print(f"❌ Error: The file '{filename}' was not found.")
    except Exception as e:
        print(f"❌ An unexpected error occurred: {e}")


if __name__ == "__main__":
    # Ensure a filename is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python replace_numbers.py <filename>")
        sys.exit(1)

    # Get the filename from the command-line arguments
    file_to_process = sys.argv[1]
    replace_numbers_in_file(file_to_process)