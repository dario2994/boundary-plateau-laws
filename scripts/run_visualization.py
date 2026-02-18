import os
import sys
import subprocess
import webbrowser
import pathlib

# Base directory of the project.
base_dir = pathlib.Path(__file__).resolve().parent.parent

# Ensure a net ID is provided
if len(sys.argv) != 2:
    sys.exit("Usage: python run_visualization.py <net_id>")

net_id = sys.argv[1]

# Validate net_id
if not net_id.isdigit():
    sys.exit("Error: net_id must be a positive integer.")

# Paths
derived_png = f"data/once_derived_images/{net_id}.png"
derived_script = "scripts/visualize_once_derived_net.py"
metrization_txt = f"data/metrized_nets/{net_id}_visualization.txt"
output_html = f"data/images3d/{net_id}.html"
spherical_script = "scripts/visualize_spherical_geodesics.py"

# Transforming the relative paths into absolute paths.
derived_png = base_dir / derived_png
derived_script = base_dir / derived_script
metrization_txt = base_dir / metrization_txt
output_html = base_dir / output_html
spherical_script = base_dir / spherical_script

# Show or generate the derived image
if not os.path.exists(derived_png):
    try:
        subprocess.run(["python", derived_script, net_id], check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error: Failed to run {derived_script}: {e}")
subprocess.run(["xdg-open", derived_png], check=False)

# Check for metrization data
if not os.path.exists(metrization_txt):
    sys.exit(f"Error: an embedding of net {net_id} was not found at {metrization_txt}")

# Run the visualization script
try:
    subprocess.run(["python", spherical_script, metrization_txt, output_html], check=True)
except subprocess.CalledProcessError as e:
    sys.exit(f"Error: Failed to run {spherical_script}: {e}")

webbrowser.open(output_html.as_uri())
