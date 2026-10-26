import argparse
import math
import pathlib
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import numpy as np


def parse_single_net_from_lines(lines, it):
    V, E, F = map(int, lines[it].split())
    it += 1
    darts = []
    for _ in range(2 * E):
        u, v = map(int, lines[it].split())
        darts.append((u, v))
        it += 1

    it += 1  # Skipping rev

    next_dart = list(map(int, lines[it].split()))
    it += 1

    cnt_concave = int(lines[it])
    it += 1

    concave_angles = list(map(int, lines[it].split())) if cnt_concave > 0 else []
    it += 1

    cnt_fake = int(lines[it])
    it += 1

    fake_vertices = list(map(int, lines[it].split())) if cnt_fake > 0 else []
    it += 1

    return (V, E, darts, next_dart, concave_angles, fake_vertices), it


def parse_one_once_derived_topology(filepath):
    with open(filepath, "r") as in_file:
        lines = [line.strip() for line in in_file if line.strip() != ""]
    net_data, _ = parse_single_net_from_lines(lines, 0)
    return net_data


def parse_many_once_derived_topologies(filepath):
    with open(filepath, "r") as in_file:
        lines = in_file.readlines()
    cnt = int(lines[0])

    once_derived_nets = []
    it = 1
    for _ in range(cnt):
        net_data, it = parse_single_net_from_lines(lines, it)
        once_derived_nets.append(net_data)
    return once_derived_nets


def draw_once_derived_net(fig, output_path, V, E, darts, next_dart, concave_angles, fake_vertices):
    # Convert the combinatorial map into a rotation system.
    rotation_system = [[] for _ in range(V)]
    darts_set = {(-1, -1)}

    for i in range(len(darts)):
        v = darts[i][0]

        # Repeated edges are not supported (the nets with repeated
        # edges are slightly changed, by another program, to make them
        # visualizable).
        assert darts[i] not in darts_set
        darts_set.add(darts[i])

        if rotation_system[v]:
            continue
        rotation_system[v].append(i)
        curr = next_dart[i]
        while curr != i:
            rotation_system[v].append(curr)
            curr = next_dart[curr]

    # Use the rotation system to define a NetworkX planar embedding.
    embedding = nx.PlanarEmbedding()
    for v in range(V):
        rs = rotation_system[v]
        embedding.add_half_edge(v, darts[rs[0]][1])

        for i in range(1, len(rs)):
            embedding.add_half_edge_ccw(v, darts[rs[i]][1], darts[rs[i - 1]][1])

    # Convert the embedding into a NetworkX graph.
    G = nx.Graph()
    for u in embedding.nodes():
        for v in embedding.adj[u]:
            G.add_edge(u, v)

    # Compute a valid planar drawing based on the embedding.
    pos = nx.combinatorial_embedding_to_pos(embedding)

    labels = {v: v for v in range(V) if v not in fake_vertices}

    nx.draw(
        G,
        pos,
        with_labels=True,
        labels=labels,
        node_color="lightblue",
        edge_color="gray",
        node_size=0,
        font_size=12,
    )

    # Highlight the angle at node b, from a->b to b->c.
    def mark_angle(a, b, c):
        if b in fake_vertices:
            return
        p1 = np.array(pos[a])
        p2 = np.array(pos[b])
        p3 = np.array(pos[c])

        v1 = p1 - p2
        v2 = p3 - p2
        angle1 = np.degrees(np.arctan2(v1[1], v1[0]))
        angle2 = np.degrees(np.arctan2(v2[1], v2[0]))
        start_angle = angle1
        angle_extent = angle2 - angle1
        if angle_extent < 0:
            angle_extent += 2 * math.pi

        arc = Arc(
            p2,
            1,
            1,
            angle=0,
            theta1=start_angle,
            theta2=start_angle + angle_extent,
            color="red",
        )
        fig.gca().add_patch(arc)

    for i in concave_angles:
        mark_angle(darts[i][0], darts[i][1], darts[next_dart[i ^ 1]][1])

    fig.gca().axis("equal")

    output_path = pathlib.Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    fig.clf()


def main():
    parser = argparse.ArgumentParser(description="Visualize once-derived topology.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--one", action="store_true", help="Visualize a single topology to an output file.")
    group.add_argument("--many", action="store_true", help="Visualize multiple topologies to an output directory.")
    parser.add_argument("--in", dest="input_path", required=True, help="Input file path.")
    parser.add_argument("--out", dest="output_path", required=True, help="Output file path (--one) or directory (--many).")

    args = parser.parse_args()

    plt.ioff()
    fig = plt.figure()

    if args.one:
        net_data = parse_one_once_derived_topology(args.input_path)
        draw_once_derived_net(fig, args.output_path, *net_data)
    else:
        nets = parse_many_once_derived_topologies(args.input_path)
        out_dir = pathlib.Path(args.output_path)
        out_dir.mkdir(parents=True, exist_ok=True)
        for idx, net_data in enumerate(nets):
            draw_once_derived_net(fig, out_dir / f"{idx}.svg", *net_data)
        print(f"Successfully generated {len(nets)} SVG vector graphics in {out_dir}")


if __name__ == "__main__":
    main()
