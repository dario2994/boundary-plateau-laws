from sympy import symbols, sin, cos, tan, Eq, pi, expand_trig, simplify, Poly, Rational, groebner, diff, resultant, sqrt, solve, factor, acos, factor_list
from sympy.parsing.sympy_parser import parse_expr
import math


sides = symbols('a b c d e')
ccos = symbols('ca cb cc cd ce')
ssin = symbols('sa sb sc sd se')
ttan = symbols('ta tb tc td te')

def ff(_x, _y):
    return cos(_x) + cos(_y) + cos(_x)*cos(_y) - 2*sin(_x)*sin(_y)

def substitute_sin_cos(eq):
    for i in range(5):
        x = sides[i]
        eq = eq.replace(cos(x), ccos[i])
        eq = eq.replace(sin(x), ssin[i])
    return eq

def cos_from_tan_half(t):
    return (1 - t**2) / (1 + t**2)

def sin_from_tan_half(t):
    return 2*t / (1 + t**2)

def substitute_tan_half(eq):
    for i in range(5):
        eq = eq.replace(cos(sides[i]), cos_from_tan_half(ttan[i]))
        eq = eq.replace(sin(sides[i]), sin_from_tan_half(ttan[i]))
    return eq

trig_system = []
for i in range(5):
    j = (i + 1) % 5
    k = (i + 3) % 5
    si, sj, sk = sides[i], sides[j], sides[k]
    trig_system.append(Eq(2*cos(sk)-Rational(1, 3), ff(si, sj)))


print("The sides a, b, c, d, e of a convex spherical pentagon with angles equal to 2pi/3 satisfy the following equations.")
print("Trigonometric System:")
for eq in trig_system:
    print(eq)
print()

poly_system = []
for eq in trig_system:
    rational_eq = substitute_tan_half(expand_trig(eq.expand())).simplify()
    rational_expr = rational_eq.lhs-rational_eq.rhs
    poly_expr = rational_expr.expand().together().as_numer_denom()[0].expand()

    poly_system.append(Poly(poly_expr))

print("The trigonometric system is equivalent to the following algebraic system. (Variable 'tV' stands for tan(V/2).)")
for eq in poly_system[0:5]:
    print(eq.as_expr(), "= 0")
print()


print("Using groebner basis in Singular (Singular < 2192.sing), this algebraic system can be shown to imply P(x, y, z) = P(z, y, x) = 0, where x=tan(a/2), y=tan(b/2), z=tan(c/2) and P is the following polynomial:")
x, y, z = symbols('x y z')
pol = 9*x**3*y**4*z**4+13*x**3*y**3*z**5-9*x**2*y**4*z**5+2*x**3*y**2*z**6-2*x**2*y**3*z**6-9*x**3*y**4*z**2-13*x**3*y**3*z**3-45*x**2*y**4*z**3+7*x**3*y**2*z**4-74*x**2*y**3*z**4+36*x*y**4*z**4+13*x**3*y*z**5-19*x**2*y**2*z**5-14*x*y**3*z**5+18*y**4*z**5+2*x**3*z**6-2*x**2*y*z**6-4*x*y**2*z**6+4*y**3*z**6-18*x**3*y**4-26*x**3*y**3*z-36*x**2*y**4*z-13*x**3*y**2*z**2-22*x**2*y**3*z**2+45*x*y**4*z**2-13*x**3*y*z**3+127*x**2*y**2*z**3-13*x*y**3*z**3+9*y**4*z**3-2*x**3*z**4+94*x**2*y*z**4-20*x*y**2*z**4-56*y**3*z**4+14*x**2*z**5-26*x*y*z**5+2*y**2*z**5-4*x*z**6+4*y*z**6-18*x**3*y**2+50*x**2*y**3+9*x*y**4-26*x**3*y*z+146*x**2*y**2*z+x*y**3*z-9*y**4*z-4*x**3*z**2+92*x**2*y*z**2-22*x*y**2*z**2-58*y**3*z**2+10*x**2*z**3-28*x*y*z**3+34*y**2*z**3+4*x*z**4+16*y*z**4-4*z**5-4*x**2*y-6*x*y**2+2*y**3-4*x**2*z-2*x*y*z+32*y**2*z+8*x*z**2+8*y*z**2-8*z**3-4*y-4*z
print(pol)
print()

print("This polynomial can be factored as the product of:")
factors = factor(pol.as_expr()).as_ordered_factors()
for f in factors:
    print("(---) ", f)

print("Observe that:")

print("The first factor is always positive.")

print("The second factor can be rewritten as:")
asym_factor = 11*x*y*z + 2*x*z**2+2*y+2*z + (9*y**2 + 2*y*z)*(x-z)
assert(factors[1] == asym_factor.expand())
print(asym_factor)

assert(factors[2].subs(x, 1).subs(z, 33) == factors[2].subs(x, 33).subs(z, 1))
print("The third factor is symmetric in x and z.")
print("Since P(x, y, z) = P(z, y, x) = 0. Either the third factor is zero or the second factor is zero and remains zero if we swap x and z. But if x >= z than the second factor is strictly positive, thus the third factor must be 0.")

Q = factors[2]
print("We have shown that Q(x, y, z) = 0, where Q is the polynomial")
print(Q)
print()

print("Now, we move our attention to the configuration appearing in net 2192. There are two pentagons sharing a side and two other sides with two quadrilaterals. Thus, if a, b, c are the three consecutive sides of one of the pentagons and b is the one shared with the other pentagon, it must hold Q(x, y, z) = Q(1/(2x), y, 1/(2z)=0.")
print("By taking the resultant of this two polynomials (with respect to y) we deduce that x, z satisfy R(x, z) = 0 where R is the polynomial:")
Q2 = (16*x**2*z**2*Q.subs(x, 1/(2*x)).subs(z, 1/(2*z))).expand()
R = resultant(Q, Q2, y)
print(R)
print()

print("The polynomial R can be factored as the product of:")

factors_R = factor(R.as_expr()).as_ordered_factors()
for f in factors_R:
    print("(---) ", f)

S = factors_R[2]
print("So, either 2xz=1, which is what we want to show, or the last factor S(x, z)=0. We know a priori that x, z are very close to 1/sqrt(2) and if we evaluate S(1/sqrt(2), 1/sqrt(2)) we obtain:")
print(S.subs(x, 1/sqrt(2)).subs(z, 1/sqrt(2)).evalf())
print("Thus, S(x, z) is not zero for the range of values we are interested in an we conclude that tan(a/2)*tan(c/2)=1/2 as desired.")