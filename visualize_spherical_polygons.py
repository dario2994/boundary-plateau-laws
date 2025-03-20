import numpy as np
import plotly.graph_objects as go
import sys
import webbrowser

def sph2cart(p):
    """Convert a 3D point on the sphere from string format to Cartesian coordinates."""
    return np.array([float(x) for x in p.strip().split()])

def great_circle_path(p, v, length, num_points=50):
    """Generates points along the geodesic starting at p, in direction v, with given arc length."""
    t = np.linspace(0, length, num_points)
    path = np.cos(t)[:, None] * p + np.sin(t)[:, None] * v
    path /= np.linalg.norm(path, axis=1)[:, None]  # Normalize to keep on sphere
    return path

def draw_spherical_polygon(data, output_file):
    """Draws the spherical polygon given a list of (p_i, v_i, len_i)."""
    fig = go.Figure()
    
    # Add sphere (Earth)
    theta = np.linspace(0, np.pi, 50)
    phi = np.linspace(0, 2 * np.pi, 50)
    x_sphere = np.outer(np.sin(theta), np.cos(phi))
    y_sphere = np.outer(np.sin(theta), np.sin(phi))
    z_sphere = np.outer(np.cos(theta), np.ones_like(phi))
    fig.add_trace(go.Surface(x=x_sphere, y=y_sphere, z=z_sphere, colorscale='Blues', opacity=0.3))
    
    # Draw geodesic edges
    for p, v, length in data:
        path = great_circle_path(p, v, length)
        fig.add_trace(go.Scatter3d(
            x=path[:, 0], y=path[:, 1], z=path[:, 2],
            mode='lines', line=dict(color='red', width=4)
        ))
    
    # Add vertices
    vertices = [p for p, _, _ in data]
    fig.add_trace(go.Scatter3d(
        x=[p[0] for p in vertices], y=[p[1] for p in vertices], z=[p[2] for p in vertices],
        mode='markers', marker=dict(size=6, color='black')
    ))
    
    # Save and show
    fig.write_html(output_file)
    webbrowser.open(output_file)

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py in.txt out.html")
        sys.exit(1)
    
    input_file, output_file = sys.argv[1], sys.argv[2]
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    n = int(lines[0].strip())
    data = []
    for i in range(1, n + 1):
        parts = lines[i].split('|')
        p = sph2cart(parts[0])
        v = sph2cart(parts[1])
        length = float(parts[2].strip())
        data.append((p, v, length))
    
    draw_spherical_polygon(data, output_file)

if __name__ == "__main__":
    main()
