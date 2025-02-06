# Function to convert the input format (adsorbate and framework) to XYZ format
def convert_to_xyz_with_framework(input_file, framework_file, output_file):
    atom_mapping = {1: 'H', 2: 'H', 3: 'O'}  # Mapping for atom types, ignoring lone pair (0)

    # Step 1: Read the water molecule file (adsorbate positions)
    xyz_lines = []
    
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Loop through each line and process adsorbate atoms
    for line in lines:
        if line.startswith('Adsorbate-atom-position'):
            parts = line.split()

            # Get the atom type and coordinates
            atom_type = int(parts[2])
            
            # Only process if the atom type is not 0 (ignore lone pairs)
            if atom_type in atom_mapping:
                x, y, z = parts[3], parts[4], parts[5]

                # Map atom type to element symbol using atom_mapping
                element = atom_mapping[atom_type]

                # Add the element and its coordinates to the XYZ format
                xyz_lines.append(f"{element} {x} {y} {z}")

    # Step 2: Read the framework atoms from the VTK file (assuming they are all carbon atoms)
    with open(framework_file, 'r') as vtk_file:
        vtk_lines = vtk_file.readlines()

    # Locate the section with POINTS in the VTK file and read the coordinates
    reading_points = False
    for line in vtk_lines:
        if line.startswith("POINTS"):
            reading_points = True
            continue  # Skip the "POINTS" line itself

        # Stop reading when encountering POINT_DATA, SCALARS, or LOOKUP_TABLE
        if line.startswith("POINT_DATA") or line.startswith("SCALARS") or line.startswith("LOOKUP_TABLE"):
            break

        # If we're in the points section, read the coordinates
        if reading_points:
            coords = line.split()
            if len(coords) == 3:  # Ensure it's a coordinate line
                x, y, z = coords
                # Add Carbon atoms (C) with the given coordinates
                xyz_lines.append(f"C {x} {y} {z}")

    # Step 3: Write to the output XYZ file
    with open(output_file, 'w') as xyz_file:
        # Write the total number of atoms (adsorbate + framework)
        xyz_file.write(f"{len(xyz_lines)}\n")
        xyz_file.write("Converted to XYZ format (including framework carbon atoms, lone pairs ignored)\n")
        
        # Write the atom positions
        for xyz_line in xyz_lines:
            xyz_file.write(xyz_line + "\n")


# Example usage:
input_file = 'restart_IC_22_AA_1.1.1_300.000000_4157'  # Input file with your adsorbate positions
framework_file = 'FrameworkAtoms.vtk'   # Framework atoms file in VTK format
output_file = 'combined_structure.xyz'  # Output XYZ file with water + framework

convert_to_xyz_with_framework(input_file, framework_file, output_file)

print(f"Conversion completed. XYZ file saved as {output_file}")

