import plotly.graph_objects as go
from math import sin, cos, sqrt

def v_add(v1, v2):
    """Computes the sum of two 3D vectors."""
    return (v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2])

def v_diff(v1, v2):
    """Computes the difference of two 3D vectors."""
    return (v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2])

def v_cross(v1, v2):
    """Computes the cross product of two 3D vectors."""
    x = v1[1] * v2[2] - v1[2] * v2[1]
    y = v1[2] * v2[0] - v1[0] * v2[2]
    z = v1[0] * v2[1] - v1[1] * v2[0]
    return (x, y, z)

def visualize_triangles(triangles, labeled_points):
    # print(triangles)
    # 3. Prepare the data for Plotly's Mesh3d object.
    #    This involves creating a list of unique vertices and mapping them to indices.
    all_points = []
    for tri in triangles:
        all_points.extend(tri)

    # Get the set of unique points
    unique_points = sorted(list(set(all_points)))
    point_to_idx = {point: i for i, point in enumerate(unique_points)}

    # Unzip unique points for the x, y, z lists
    x_coords = [p[0] for p in unique_points]
    y_coords = [p[1] for p in unique_points]
    z_coords = [p[2] for p in unique_points]

    # Create the i, j, k index lists for the triangles
    i_indices, j_indices, k_indices = [], [], []
    for tri in triangles:
        i_indices.append(point_to_idx[tri[0]])
        j_indices.append(point_to_idx[tri[1]])
        k_indices.append(point_to_idx[tri[2]])

    # 4. Create the figure
    fig = go.Figure(data=[
        go.Mesh3d(
            x=x_coords,
            y=y_coords,
            z=z_coords,
            i=i_indices,
            j=j_indices,
            k=k_indices,
            opacity=0.8,
            color='cyan',
            flatshading=True # Gives a faceted look
        )
    ])

    # Improve the layout
    fig.update_layout(
        title_text="Visualization of 3D Triangles",
        scene=dict(
            xaxis=dict(title='X-axis'),
            yaxis=dict(title='Y-axis'),
            zaxis=dict(title='Z-axis'),
            aspectmode='data' # This ensures the aspect ratio is not distorted
        )
    )

    edge_x, edge_y, edge_z = [], [], []
    for tri in triangles:
        # Add edges p1->p2, p2->p3, p3->p1
        edge_x.extend([tri[0][0], tri[1][0], tri[2][0], tri[0][0], None])
        edge_y.extend([tri[0][1], tri[1][1], tri[2][1], tri[0][1], None])
        edge_z.extend([tri[0][2], tri[1][2], tri[2][2], tri[0][2], None])

    # Edges
    edge_trace = go.Scatter3d(
        x=edge_x,
        y=edge_y,
        z=edge_z,
        mode='lines',
        line=dict(color='black', width=5),
        name='Edges'
    )
    fig.add_trace(edge_trace)

    # Points
    if labeled_points:
        point_x = [p[0][0] for p in labeled_points]
        point_y = [p[0][1] for p in labeled_points]
        point_z = [p[0][2] for p in labeled_points]
        point_labels = [p[1] for p in labeled_points]

        points_trace = go.Scatter3d(
            x=point_x,
            y=point_y,
            z=point_z,
            mode='markers+text',
            text=point_labels,
            marker=dict(
                color='red',
                size=5
            ),
            textposition='top center',
            name='Labeled Points'
        )
        fig.add_trace(points_trace)



    fig.show()


def triangle_area_squared(p, q, r):
    pq = (q[0] - p[0], q[1] - p[1], q[2] - p[2])
    pr = (r[0] - p[0], r[1] - p[1], r[2] - p[2])

    dot_pq_pq = pq[0]**2 + pq[1]**2 + pq[2]**2
    dot_pr_pr = pr[0]**2 + pr[1]**2 + pr[2]**2
    dot_pq_pr = pq[0] * pr[0] + pq[1] * pr[1] + pq[2] * pr[2]

    return (dot_pq_pq * dot_pr_pr - dot_pq_pr**2) / 4

def triangle_area(p, q, r):
    return sqrt(triangle_area_squared(p, q, r))

def skew_quadrilateral_area1(a, b, c, d):
    return min(triangle_area(a, b, c) + triangle_area(a, c, d),
               triangle_area(a, b, d) + triangle_area(b, c, d))


import numpy as np
from scipy.integrate import dblquad

def compute_surface_area(v1, v2, v3):
    """
    Computes the area of a bilinear patch given its defining vectors.

    The area is calculated by the integral over the unit square:
    Area = ∫[0,1] ∫[0,1] ||v1 + s*v2 + t*v3|| ds dt

    Args:
        v1 (tuple or np.ndarray): The constant vector component.
        v2 (tuple or np.ndarray): The vector component scaled by 's'.
        v3 (tuple or np.ndarray): The vector component scaled by 't'.

    Returns:
        float: The calculated area of the surface patch.
    """
    # Convert inputs to numpy arrays for vector operations
    vec1 = np.asarray(v1)
    vec2 = np.asarray(v2)
    vec3 = np.asarray(v3)

    # Define the integrand function locally
    # Note: dblquad requires the function signature to be f(y, x),
    # so we use f(t, s).
    def integrand(t, s):
        v = vec1 + s * vec2 + t * vec3
        return np.linalg.norm(v)

    # Perform the numerical integration
    # The result is a tuple (area, estimated_error)
    area, _ = dblquad(
        integrand, # function
        0,         # s-limit from
        1,         # s-limit to
        0,         # t-limit from
        1          # t-limit to
    )

    return area

def skew_quadrilateral_area2(a, b, c, d):
    v0 = v_diff(v_add(a, c), v_add(b, d))
    v1 = v_cross(v_diff(b, a), v_diff(d, a))
    v2 = v_cross(v_diff(b, a), v0)
    v3 = v_cross(v0, v_diff(d, a))

    return compute_surface_area(v1, v2, v3)


def skew_quadrilateral_diff(a, b, c, d):
    return skew_quadrilateral_area1(a, b, c, d) - skew_quadrilateral_area2(a, b, c, d)

def sym_x(x, y, z):
    return (-x, y, z)

def sym_z(x, y, z):
    return (x, y, -z)

def make_computation(s, t, q):
    def rotate(x, y, z):
        return (x*cos(q)-y*sin(q), y*cos(q)+x*sin(q), z)

    pin0 = (0, 0, sqrt(3))
    pin1 = sym_z(*pin0)
    A = (0, 0, -s)
    C = sym_z(*A)
    B = (0, 0, 0)
    D = (0, t, 0)

    P = [(0, 0, 0) for _ in range(8)]
    # Face of the cube
    P[1] = (1, sqrt(2), 0)
    P[4] = (1, 0, sqrt(2))
    P[6] = sym_z(*P[4])
    P[0] = (1, -sqrt(2), 0)

    for i in [1, 4, 6, 0]:
        P[i] = rotate(*P[i])

    for (i, j) in [(4, 5), (0, 3), (6, 7), (1, 2)]:
        P[j] = sym_x(*P[i])

    triangles = [
        (C, D, A),
        (D, P[1], P[2]),
        (B, P[0], P[3]),
        (C, P[4], pin0),
        (C, P[5], pin0),
        (A, P[6], pin1),
        (A, P[7], pin1),
        (C, D, P[1]),
        (C, D, P[2]),
        (C, P[4], P[1]),
        (C, P[5], P[2]),
        (C, B, P[0]),
        (C, B, P[3]),
        (P[4], C, P[0]),
        (P[5], C, P[3]),
        (B, A, P[0]),
        (B, A, P[3]),
        (A, P[0], P[6]),
        (A, P[3], P[7]),
        (A, P[1], D),
        (A, P[2], D),
        (P[1], A, P[6]),
        (P[2], A, P[7]),
    ]

    visualize_triangles(triangles,
        [(P[i], str(i)) for i in range(len(P))] +
        [(pp[0], pp[1]) for pp in [(A, 'A'), (B, 'B'), (C, 'C'), (D, 'D'), (pin0, 'pin0'), (pin1, 'pin1')]])

    # diff1 = skew_quadrilateral_diff(C, D, P[1], P[4])
    # diff2 = skew_quadrilateral_diff(D, A, P[6], P[1])
    # diff3 = skew_quadrilateral_diff(C, B, P[0], P[4])
    # diff4 = skew_quadrilateral_diff(A, B, P[0], P[6])
    # gain = 2 * (diff1 + diff2 + diff3 + diff4)
    # print(diff1, diff2, diff3, diff4, gain)

    tot_area = 0
    for triangle in triangles[0:7]:
        area = triangle_area(*triangle)
        tot_area += area

    quad1 = skew_quadrilateral_area2(C, D, P[1], P[4])
    quad2 = skew_quadrilateral_area2(D, A, P[6], P[1])
    quad3 = skew_quadrilateral_area2(C, B, P[0], P[4])
    quad4 = skew_quadrilateral_area2(A, B, P[0], P[6])
    tot_area += 2 * (quad1 + quad2 + quad3 + quad4)

    return tot_area

q_val = 0.1
s_val = 0.1
t_val = 0.15

area_cone = make_computation(s_val, t_val, q_val)
print("Cone:", area_cone)
exit(0)

print("Area competitors:")
for s_val in np.arange(0.0, 1, 0.02):
    t_best = -1
    lowest_area = 1e9
    for t_val in np.arange(0.0, 2, 0.02):
        area_competitor_val = make_computation(s_val, t_val, q_val)
        if area_competitor_val < lowest_area:
            t_best = t_val
            lowest_area = area_competitor_val
    print(s_val, t_best, lowest_area - area_cone)


# print("Competitor:", make_computation(s_val, t_val, q_val))

# visualize_triangles(evaluate_triangles(s_val, t_val, q_val),
#     [(evaluate_point(s_val, t_val, P[i]), str(i)) for i in range(len(P))] +
#     [(evaluate_point(s_val, t_val, pp[0]), pp[1]) for pp in [(A, 'A'), (B, 'B'), (C, 'C'), (D, 'D')]] +
#     [(evaluate_point(s_val, t_val, pin0), 'pin0'),
#      (evaluate_point(s_val, t_val, pin1), 'pin1')])

# def print_comparison(s_val, t_val):
#     P_val = [evaluate_point(s_val, t_val, p) for p in P]
#     A_val, B_val, C_val, D_val, pin0_val, pin1_val = \
#         [evaluate_point(s_val, t_val, p) for p in [A, B, C, D, pin0, pin1]]
#     print(skew_quadrilateral_area1(C_val, D_val, P_val[1], P_val[4]),
#           skew_quadrilateral_area2(C_val, D_val, P_val[1], P_val[4]))


exit(0)


# tot_area = expand_trig(tot_area.expand().simplify())

# area_cone = tot_area.subs(s, 0).subs(t, 0)
# area_cone = expand_trig(area_cone.expand().simplify())

# area_competitor = expand_trig(tot_area.expand().simplify()).simplify()

# area_cone = area_cone.subs(q, q_val)
# area_competitor = area_competitor.subs(q, q_val)

# print("Area cone:")
# print(area_cone.evalf())

# print(area_competitor)

# def fast_area_competitor(_s, _t):
#     return _s*_t + 2*math.sqrt(3)*_s + 4*math.sqrt(2)*math.sqrt(_s**2 - 2*math.sqrt(2)*_s + 4) + 2*math.sqrt(_s**2 - 2*math.sqrt(3)*_s + 3) + math.sqrt(-0.399333666587313*math.sqrt(2)*_t**2 + 4.00999166944395*_t**2 - 7.97004371824188*math.sqrt(2)*_t + 1.19450630911192*_t - 0.397338661590122*math.sqrt(2) + 7.93023302244434)/2 + 2*math.sqrt(_s**2*_t**2 - 1.99750052078993*math.sqrt(2)*_s**2*_t - 0.0999583385413567*_s**2*_t + 3*_s**2 - 0.0998334166468282*math.sqrt(2)*_t**2 + 1.00249791736099*_t**2) + math.sqrt(0.397338661590122*math.sqrt(2) + 7.93023302244434)/2

# # print(area_competitor.subs(s, s_val).subs(t, t_val).evalf(), fast_area_competitor(s_val, t_val))


# print("Area competitors:")
# for s_val in np.arange(0.0, 0.5, 0.01):
#     t_best = -1
#     lowest_area = 1e9
#     for t_val in numpy.arange(0.0, 1, 0.01):
#         area_competitor_val = fast_area_competitor(s_val, t_val)
#         if area_competitor_val < lowest_area:
#             t_best = t_val
#             lowest_area = area_competitor_val
#     print(s_val, t_best, lowest_area)
# exit(0)
# import numpy

# for s_val in numpy.arange(0.001, 0.3, 0.005):
#     for t_val in numpy.arange(0.001, 0.3, 0.005):
#         area_competitor_val = area_competitor.subs(s, s_val).subs(t, t_val).subs(q, q_val).evalf()
#         # print(area_val)
#         # print(s_val, t_val)
#         if area_competitor_val < area_cone_val:
#             print("Area competitor:")
#             print(s_val, t_val, area_competitor_val)
