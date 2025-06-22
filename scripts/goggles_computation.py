from sympy import symbols, sin, cos, tan, Eq, pi, expand_trig, simplify, Poly, Rational, groebner, diff, resultant, sqrt
from sympy.parsing.sympy_parser import parse_expr


a, b, c, alpha, beta = symbols('a b c alpha beta')
tb, tc, talpha, tbeta = symbols('tb tc talpha tbeta')

def cos_from_tan_half(t):
    return (1 - t**2) / (1 + t**2)

def sin_from_tan_half(t):
    return 2*t / (1 + t**2)

def substitute_tan_half(eq):
    for (x, tx) in [(a, ta), (b, tb), (c, tc), (alpha, talpha), (beta, tbeta)]:
        eq = eq.replace(cos(x), cos_from_tan_half(tx))
        eq = eq.replace(sin(x), sin_from_tan_half(tx))
    return eq



beta_prime = beta-pi/3     # Additional equation 1.
ta = 1/(2*tb)              # Additional equation 2.
trig_system = [
    Eq(cos(c), cos(a) * cos(b) + cos(2*pi/3)*sin(a)*sin(b)),
    Eq(sin(beta), sin(b)*sin(2*pi/3)/sin(c)),
    Eq(cos(alpha), cos(c)*sin(beta_prime)*sin(pi-beta_prime)-cos(beta_prime)*cos(pi-beta_prime)),
]

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

print("Polynomial System:")
for eq in poly_system:
    print(eq)
    print()

# print(poly_system[1].reduce(poly_system[1]))
# print()
#############################################################
# The Groebner basis calculation requires a few minutes.
# G = groebner(poly_system, tc, tbeta, talpha, tb, order='lex')
# pols_in_tb_talpha = [g for g in G if g.free_symbols <= {tb, talpha}]
# assert(len(pols_in_tb_talpha) == 1)
# pol = Poly(pols_in_tb_talpha[0].as_expr())
#############################################################
pol = Poly(parse_expr("talpha**4*tb**14 + 9/2*talpha**4*tb**12 + 15/2*talpha**4*tb**10 + 49/8*talpha**4*tb**8 + 45/16*talpha**4*tb**6 + 3/4*talpha**4*tb**4 + 1/16*talpha**4*tb**2 - 2/3*talpha**2*tb**14 - 7/3*talpha**2*tb**12 - 4*talpha**2*tb**10 - 13/3*talpha**2*tb**8 - 35/12*talpha**2*tb**6 - 19/16*talpha**2*tb**4 - 7/24*talpha**2*tb**2 - 1/48*talpha**2 + 1/9*tb**14 + 5/18*tb**12 + 1/6*tb**10 - 5/72*tb**8 - 13/144*tb**6 - 1/48*tb**4 + 1/144*tb**2 + 1/144"))
#############################################################
der_pol = diff(pol, tb)

print("Equation relating tan(b/2) and tan(alpha/2):")
print(pol)
print()

print("Derivative of the previous equation with respect to tb:")
print(der_pol)
print()


res_tb = Poly(resultant(pol, der_pol, talpha))
print("Any value of tan(b/2) such at this value tan(alpha/2) has null derivative is a root of the following polynomial")
print(res_tb)
print()

tb_values = res_tb.real_roots()
print("The possible values of tan(b/2) are:")
print(res_tb.real_roots())
