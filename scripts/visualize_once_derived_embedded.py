import argparse
import pathlib
import numpy as np
import plotly.graph_objects as go


def great_circle_path(p, v, length, num_points=50):
    """Generates points along the geodesic starting at p, in direction v, with given arc length."""
    t = np.linspace(0, length, num_points)
    path = np.cos(t)[:, None] * p + np.sin(t)[:, None] * v
    path /= np.linalg.norm(path, axis=1)[:, None]
    return path


def parse_metrization_file(filepath):
    """Parses a metrized net file into vertices and edges."""
    with open(filepath, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    v, e = map(int, lines[0].split())
    vertices = []
    edges = []

    for i in range(1, v + 1):
        values = lines[i].split()
        assert len(values) == 4
        x, y, z = float(values[0]), float(values[1]), float(values[2])
        label = values[3]
        vertices.append(([x, y, z], label))

    for i in range(v + 1, e + v + 1):
        values = list(map(float, lines[i].split()))
        assert len(values) == 7
        p = np.array(values[0:3])
        vec = np.array(values[3:6])
        length = values[6]
        edges.append((p, vec, length))

    return vertices, edges


def build_3d_net_traces(vertices, edges):
    """Builds serialized Plotly-compatible traces for the spherical net."""
    theta = np.linspace(0, np.pi, 40)
    phi = np.linspace(0, 2 * np.pi, 40)
    x_sphere = (np.outer(np.sin(theta), np.cos(phi))).tolist()
    y_sphere = (np.outer(np.sin(theta), np.sin(phi))).tolist()
    z_sphere = (np.outer(np.cos(theta), np.ones_like(phi))).tolist()

    traces = [
        {
            "type": "surface",
            "x": x_sphere,
            "y": y_sphere,
            "z": z_sphere,
            "colorscale": "Blues",
            "opacity": 0.25,
            "showscale": False,
            "hoverinfo": "skip",
        }
    ]

    # Single combined trace for all geodesic edges
    edge_x, edge_y, edge_z = [], [], []
    for p, v, length in edges:
        path = great_circle_path(p, v, length)
        edge_x.extend(path[:, 0].tolist() + [None])
        edge_y.extend(path[:, 1].tolist() + [None])
        edge_z.extend(path[:, 2].tolist() + [None])

    traces.append(
        {
            "type": "scatter3d",
            "mode": "lines",
            "x": edge_x,
            "y": edge_y,
            "z": edge_z,
            "line": {"color": "red", "width": 4},
            "hoverinfo": "skip",
            "showlegend": False,
        }
    )

    # Vertices
    traces.append(
        {
            "type": "scatter3d",
            "mode": "markers+text",
            "x": [pos[0] for pos, _ in vertices],
            "y": [pos[1] for pos, _ in vertices],
            "z": [pos[2] for pos, _ in vertices],
            "text": [label for _, label in vertices],
            "textposition": "top center",
            "marker": {"size": 4, "color": "green"},
            "textfont": {"color": "black", "size": 13},
            "hoverinfo": "skip",
            "showlegend": False,
        }
    )

    return traces


def draw_geodesic_net(vertices, edges, output_file):
    """Renders traces to a standalone HTML file."""
    traces = build_3d_net_traces(vertices, edges)
    fig = go.Figure(data=traces)
    fig.update_layout(
        hovermode=False,
        scene=dict(
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title="", showspikes=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title="", showspikes=False),
            zaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title="", showspikes=False),
            aspectmode="cube",
        ),
        margin=dict(l=0, r=0, b=0, t=0),
    )
    output_path = pathlib.Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_path))


def main():
    parser = argparse.ArgumentParser(description="Visualize embedded once-derived net.")
    parser.add_argument("--in", dest="input_file", required=True, help="Input metrization .txt file.")
    parser.add_argument("--out", dest="output_file", required=True, help="Output standalone .html file.")
    args = parser.parse_args()

    vertices, edges = parse_metrization_file(args.input_file)
    draw_geodesic_net(vertices, edges, args.output_file)


if __name__ == "__main__":
    main()
