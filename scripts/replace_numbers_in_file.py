import sys
import re

def replace_numbers_in_file(filename):
    """
    Replaces specific numbers in a file according to a predefined mapping.
    The replacement happens only for whole numbers.
    """
    # Define the mapping from old numbers to new numbers
    replacements = {
        '28': '28', '100': '99', '160': '158', '184': '182', '252': '250', '322': '319', '328': '325', '350': '346', '365': '361', '426': '421', '474': '468', '490': '484', '526': '520', '621': '614', '690': '681', '722': '712', '747': '736', '790': '779', '839': '827', '924': '909', '935': '920', '954': '938', '963': '-36178914', '993': '-36178914', '1048': '1029', '1118': '1096', '1122': '1100', '1147': '1122', '1160': '-36178914', '100': '99', '101': '100', '11': '11', '39': '39', '722': '712', '360': '356', '45': '45', '1122': '1100', '1067': '1047'
    }

    legacy_replacements = {'72': '28', '278': '99', '434': '158', '521': '182', '691': '250', '936': '319', '949': '325', '1011': '346', '1034': '361', '1223': '421', '1351': '468', '1405': '484', '1514': '520', '1829': '614', '2017': '681', '2119': '712', '2192': '736', '2312': '779', '2454': '827', '2676': '909', '2700': '920', '2769': '938', '2799': '-36178914', '2891': '-36178914', '3048': '1029', '3222': '1096', '3233': '1100', '3295': '1122', '3338': '-36178914'}

    # Update legacy replacements.
    # legacy_replacements = {key: replacements[value] for key, value in legacy_replacements.items()}
    # print(legacy_replacements)

    try:
        # Read the original content of the file
        with open(filename, 'r') as file:
            content = file.read()

        # Keep track of original content to see if changes were made
        original_content = content

        # Iterate through the replacement map and apply changes
        for old_num, new_num in legacy_replacements.items():
            # Use regex with word boundaries (\b) to match whole numbers only.
            # This ensures '72' is matched, but '172' or '72a' are not.
            old_pattern = r'\b' + old_num + r'\b'
            # new_pattern = new_num
            new_pattern = "{\\\color{red} " + new_num + " Warning}"
            content = re.sub(old_pattern, new_pattern, content)

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