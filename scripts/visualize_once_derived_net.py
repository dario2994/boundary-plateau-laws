import gc
import math
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import numpy as np
import sys
import pathlib
import os

# Base directory of the project.
base_dir = pathlib.Path(__file__).resolve().parent.parent
once_derived_for_visualization = base_dir / 'data/once_derived_for_visualization.txt'
once_derived_images_dir = base_dir / 'data/once_derived_images/'

# Ensure a net ID is provided
if len(sys.argv) != 2:
    sys.exit("Usage: python visualize_once_derived_net.py <net_id>")

net_id = sys.argv[1]

# Validate net_id
if not net_id.isdigit():
    sys.exit("Error: net_id must be a positive integer.")

def parse_once_derived_nets():
    lines = []
    with open(once_derived_for_visualization) as in_file:
        lines = in_file.readlines()
    cnt = int(lines[0])

    once_derived_nets = []
    it = 1
    for z in range(cnt):
        V, E, F = map(int, lines[it].split())
        it += 1
        darts = []
        for i in range(2*E):
            u, v = map(int, lines[it].split())
            darts.append((u, v))
            it += 1

        it += 1 # Skipping rev

        next_dart = list(map(int, lines[it].split()))
        it += 1

        cnt_concave = int(lines[it])
        it += 1

        concave_angles = list(map(int, lines[it].split()))
        it += 1

        cnt_fake = int(lines[it])
        it += 1
        fake_vertices = list(map(int, lines[it].split()))
        it += 1

        once_derived_nets.append((V, E, darts, next_dart, concave_angles, fake_vertices))
    return once_derived_nets


def draw_once_derived_net(fig, _id, V, E, darts, next_dart, concave_angles, fake_vertices):
    print('Drawing once-derived net ' + str(_id) +'.')

    # Convert the combinatorial map in a rotation system.
    rotation_system = [[] for _ in range(V)]
    darts_set = {(-1, -1)}

    for i in range(len(darts)):
        v = darts[i][0]

        # Repeated edges are not supported (the nets with repeated
        # edges are slightly changed, by another program, to make them
        # visualizable).
        assert(darts[i] not in darts_set)
        darts_set.add(darts[i])

        if rotation_system[v]:
            continue
        rotation_system[v].append(i)
        curr = next_dart[i]
        while curr != i:
            rotation_system[v].append(curr)
            curr = next_dart[curr]

    # Use the rotation system to define a networkx planar embedding.
    embedding = nx.PlanarEmbedding()
    for v in range(V):
        rs = rotation_system[v]
        embedding.add_half_edge(v, darts[rs[0]][1])

        for i in range(1, len(rs)):
            embedding.add_half_edge_ccw(v, darts[rs[i]][1], darts[rs[i-1]][1])

    # Convert the embedding into a NetworkX graph
    G = nx.Graph()
    for u in embedding.nodes():
        for v in embedding.adj[u]:
            G.add_edge(u, v)

    # Compute a valid planar drawing based on the embedding
    pos = nx.combinatorial_embedding_to_pos(embedding)

    labels = {}
    for v in range(V):
        if v not in fake_vertices:
            labels[v] = v

    nx.draw(G, pos, with_labels=True, labels=labels, node_color='lightblue', edge_color='gray', node_size=0, font_size=12)


    # Highlight the angle at node a, from a->b to b->c
    def mark_angle(a, b, c):
        if b in fake_vertices:
            return
        node = b
        # Positions of the three nodes.
        p1 = np.array(pos[a])
        p2 = np.array(pos[b])
        p3 = np.array(pos[c])


        v1 = p1 - p2
        v2 = p3 - p2
        angle1 = np.degrees(np.arctan2(v1[1], v1[0]))  # Angle of 1->2
        angle2 = np.degrees(np.arctan2(v2[1], v2[0]))  # Angle of 2->3
        start_angle = angle1
        angle_extent = angle2 - angle1
        if angle_extent < 0:
            angle_extent += 2*math.pi

        # Draw an arc at node 2 to represent the angle
        arc = Arc(p2, 1, 1, angle=0, theta1=start_angle, theta2=start_angle + angle_extent, color='red')
        fig.gca().add_patch(arc)

    for i in concave_angles:
        mark_angle(darts[i][0], darts[i][1], darts[next_dart[i^1]][1])

    fig.gca().axis('equal')
    fig.gca().set_title("once-derived net: " + str(_id))

    # ~ plt.show()
    os.makedirs(once_derived_images_dir, exist_ok=True)
    fig.savefig(once_derived_images_dir / f'{_id}.png', dpi=300, bbox_inches='tight')
    fig.clf()


once_derived_nets = parse_once_derived_nets()

# ~ survivors = []
# ~ with open('data/once_derived_survivors.txt') as in_file:
    # ~ survivors = list(map(int, in_file.readlines()[1].split()))

plt.ioff()
# We reuse the same figure as a workaround for a memory leak in
# matplotlib.
fig = plt.figure()

draw_once_derived_net(fig, net_id, *(once_derived_nets[int(net_id)]))
