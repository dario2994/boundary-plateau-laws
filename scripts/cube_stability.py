# TODEL: This file is most likely broken/wrong. I don't want to delete it as
#        it could be a good starting point of a working version.

import sympy.vector
from sympy import sin, cos, sqrt, symbols, collect

N = sympy.vector.CoordSys3D('N')

def vector3d(x, y, z):
    return x * N.i + y * N.j + z * N.k

s, t = symbols('s t')
l = symbols('l') # Side of the small square
q = symbols('q') # Rotation angle

# expr = (8*t**2 - 16*t + 8)*(-2*sqrt(2)*l*s*sin(q)/(4*t - 4) - 2*l*s*cos(q)/(4*t - 4) + 2*l*s/(4*t - 4) - sqrt(2)*l*t*sin(q)/(4*t - 4) - 3*l*t*cos(q)/(4*t - 4) - 3*l*t/(4*t - 4) + 2*sqrt(2)*l*sin(q)/(4*t - 4) + 2*l*cos(q)/(4*t - 4))
# print(expr.cancel())

# exit(0)

# Taylor expand with respect to l
def quadrilateral_area(A, B, C, D):
    v1 = (B-A)^(D-A)
    v2 = (B-A)^(A-B-D+C)
    v3 = (A-B-D+C)^(D-A)
    vec = v1 + s * v2 + t * v3
    # We divide and multiply by 1-t because sympy would not be able to simplify
    # sqrt((1-t)²) because it does not know the sign of 1-t.
    vec /= (1-t)
    f = sqrt(vec&vec).series(l, 0, 3).removeO()*(1-t)
    print(f.simplify().cancel().simplify())

    # print(f)

def Quad():
    A = vector3d(0, 1, 1/sqrt(2))
    B = vector3d(1, 0, 1/sqrt(2))
    C = vector3d(1, 0, 0)
    D = vector3d(0, cos(q), sin(q))

    return quadrilateral_area(A, B, C*l, D*l)

def triangle_area(A, B, C):
    return ((B-A)^(C-A))/2

def Tri():
    A = vector3d(0, 1, 1/sqrt(2))
    B = vector3d(0, cos(q), sin(q))
    C = vector3d(0, 1, -1/sqrt(2))

    return triangle_area(A, B*l, C)

Quad()