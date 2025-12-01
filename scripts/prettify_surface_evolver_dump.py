import re
import sys

def parse_evolver_dmp(input_filepath, output_filepath):
    """
    Parses a Surface Evolver .dmp file to extract the positions of vertices
    and the simplicial structure (edges and faces).

    Args:
        input_filepath (str): The path to the input .dmp file.
        output_filepath (str): The path where the formatted output will be saved.
    """
    try:
        with open(input_filepath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: Input file not found at '{input_filepath}'")
        return

    vertices = []
    edges = []
    faces = []

    # Mappings from 1-based file IDs to 0-based list indices
    vertex_map = {}
    edge_map = {}

    current_section = None

    for line in lines:
        stripped = line.strip().lower()
        # print(stripped)

        # Skip empty lines and comments
        if not stripped or stripped.startswith('//'):
            continue

        if stripped.startswith("vertices "):
            current_section = "vertices"
            continue
        elif stripped == "edges":
            current_section = "edges"
            continue
        elif stripped.startswith("faces "):
            current_section = "faces"
            continue
        elif stripped.startswith("bodies") or  stripped.startswith("read"):
            current_section = None
            continue

        # Regex to find numbers, including scientific notation
        parts = re.findall(r'-?[\d\.]+(?:e-?\d+)?', stripped)
        if current_section == "vertices":
            # Boundary vertex
            if 'boundary' in stripped and '/*' in stripped:
                v_id = int(parts[0])
                assert(len(parts) == 6)
                coords = [float(p) for p in parts[3:6]]
                vertex_map[v_id] = len(vertices)
                vertices.append(coords)
            # Normal vertex
            else:
                # if len(parts) != 4:
                #     print(stripped, parts, 'original' in stripped)
                assert(len(parts) == 4 or (len(parts) == 5 and 'original' in stripped))
                v_id = int(parts[0])
                coords = [float(p) for p in parts[1:4]]
                vertex_map[v_id] = len(vertices)
                vertices.append(coords)

        elif current_section == "edges":
            assert(len(parts) in [3, 4, 5])
            e_id, v1_id, v2_id = map(int, parts[:3])
            assert(v1_id != v2_id)
            edge_map[e_id] = len(edges)
            edges.append((vertex_map[v1_id], vertex_map[v2_id]))

        elif current_section == "faces":
            assert(len(parts) in [5, 6])
            triangle_vertices = set()
            for e_str in parts[1:4]:
                e_id = edge_map[abs(int(e_str))]
                triangle_vertices.add(edges[e_id][0])
                triangle_vertices.add(edges[e_id][1])
            assert(len(triangle_vertices) == 3)
            faces.append(tuple(triangle_vertices))

    # Write the formatted data to the output file
    try:
        with open(output_filepath, 'w') as f:
            # Header line
            f.write(f"{len(vertices)} {len(faces)}\n")

            # Vertex coordinates
            for v_coords in vertices:
                f.write(f"{v_coords[0]} {v_coords[1]} {v_coords[2]}\n")

            # Faces
            for e1, e2, e3 in faces:
                f.write(f"3 {e1} {e2} {e3}\n")

        print(f"Successfully parsed file and wrote data to '{output_filepath}'.")
        print(f"Summary: {len(vertices)} vertices, {len(edges)} edges, {len(faces)} faces.")

    except IOError:
        print(f"Error: Could not write to output file at '{output_filepath}'")



if __name__ == '__main__':
    import sys

    if len(sys.argv) != 3:
        print("Usage: python prettify_surface_evolver_dump.py <input_filename> <output_filename>")
    else:
        input_filename = sys.argv[1]
        output_filename = sys.argv[2]
        parse_evolver_dmp(input_filename, output_filename)