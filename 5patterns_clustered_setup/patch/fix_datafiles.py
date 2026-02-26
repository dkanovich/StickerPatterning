import os
import sys

def fix(data):
    # Read the file content
    with open(data, 'r') as file:
        lines = file.readlines()

    modified_lines = []
    for line in lines:
        # Check for coordinate lines (e.g., xlo xhi, ylo yhi, zlo zhi)
        if any(coord in line for coord in ['xlo xhi', 'ylo yhi', 'zlo zhi']):
            parts = line.split()
            if len(parts) == 4:  # Ensure the line has the correct format
                try:
                    # Increase the coordinates by 100
                    parts[0] = str(float(parts[0]) - 100)
                    parts[1] = str(float(parts[1]) + 100)
                    modified_line = ' '.join(parts) + '\n'
                    modified_lines.append(modified_line)
                    continue
                except ValueError:
                    pass  # If conversion fails, keep the line unchanged

        # Add other modifications if needed
        elif '2  bond types' in line:
            print(f"fixed data file {data}")
            modified_line = line.replace('2', '3')  # Change '2' to '3'
            modified_lines.append(modified_line)
            modified_lines.append('50 extra bond per atom\n')  # Insert new line
            continue

        # If no modifications are needed, keep the line unchanged
        modified_lines.append(line)

    # Overwrite the file with the modified content
    with open(data, 'w') as file:
        file.writelines(modified_lines)

if len(sys.argv) != 2:
    print("Usage: python script.py <datafile>")
    sys.exit(1)

datafile = sys.argv[1]
fix(datafile)
print("Fixed data file")