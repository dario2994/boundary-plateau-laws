from sympy import symbols, sin, cos, tan, Eq, pi, expand_trig, simplify, Poly, Rational, groebner, diff, resultant, sqrt, solve, factor, acos, factor_list
from sympy.parsing.sympy_parser import parse_expr

def sin_rule(A, B, a, b):
    return Eq(sin(A)*sin(b), sin(B)*sin(a))

def cos_rule_angles(A, B, C, a):
    return Eq(cos(A), -cos(B)*cos(C)+sin(B)*sin(C)*cos(a))

def cos_rule_sides(a, b, c, A):
    return Eq(cos(a), cos(b)*cos(c) + sin(b)*sin(c)*cos(A))

def cot_rule(A, c, B, a):
    return Eq(cos(c)*cos(B)*sin(a)*sin(A), sin(A)*cos(a)*sin(c)-sin(a)*cos(A)*sin(B))


BAC, AC, AD, ADC, EX, DEX = symbols('BAC AC AD ADC EX DEX')
ccos = symbols('c0 c1 c2 c3 c4 c5')
ssin = symbols('s0 s1 s2 s3 s4 s5')

ellT = acos(-Rational(1, 3))
ellT2 = acos(1/sqrt(3)) # = ellT/2
BAC_val, AC_val, AD_val, ADC_val, EX_val, DEX_val = pi/2, ellT2, ellT2, 2*pi/3, ellT, pi/3


variables = [BAC, AC, AD, ADC, EX, DEX]
var_vals = [BAC_val, AC_val, AD_val, ADC_val, EX_val, DEX_val]

var_num = len(variables)

def substitute_sin_cos(eq):
    for i in range(var_num):
        x = variables[i]
        eq = eq.replace(cos(x), ccos[i])
        eq = eq.replace(sin(x), ssin[i])
    return eq

def substitute_values(eq):
    for i in range(var_num):
        eq = eq.subs(variables[i], var_vals[i])
    return expand_trig(eq.expand()).simplify()

trig_system = [[], []]

# Four equations, 4 variables
trig_system[0].append(cos_rule_angles(pi/3, pi/3, BAC, AC))
trig_system[0].append(sin_rule(ADC, pi/3, AC, AD))
trig_system[0].append(cos_rule_angles(pi/3, pi/2-BAC, ADC, AD))
trig_system[0].append(Eq(2*cos(AD)*sin(ADC), 1)) # This additional equation is explained later on.

# Three equations, 4 variables
trig_system[1].append(sin_rule(ADC, pi/3, EX, pi-2*AD))
trig_system[1].append(cot_rule(ADC, pi-2*AD, DEX, EX))
trig_system[1].append(cot_rule(pi/3, EX, pi/2 - DEX, pi/2))


for i in range(2):
    print("Trigonometric System", i, ":")
    for eq in trig_system[i][0:3]:
        print(eq)
        assert(substitute_values(eq))
        print()

sqrt3 = symbols('sqrt3')
poly_system = [[], []]
for i in range(2):
    for eq in trig_system[i]:
        rational_eq = substitute_sin_cos(expand_trig(eq.expand())).simplify()
        rational_eq = rational_eq.subs(sqrt(3), sqrt3)
        rational_expr = rational_eq.lhs-rational_eq.rhs
        poly_expr = rational_expr.expand().together().as_numer_denom()[0].expand()

        poly_system[i].append(Poly(poly_expr))
    poly_system[i].append(Poly(sqrt3**2-3))


for i in range(0, 4):
    poly_system[0].append(Poly(ssin[i]**2+ccos[i]**2-1))
for i in range(2, 6):
    poly_system[1].append(Poly(ssin[i]**2+ccos[i]**2-1))

print("The two trigonometric systems are equivalent, to the following algebraic systems. (Variables starting with s and c are respectively sines and cosines of the variables appearing above.)")
for i in range(2):
    print("Polynomial System", i, ":")
    for eq in poly_system[i][0:3]:
        print(eq.as_expr(), "= 0")
    print()

# print(cos(var_vals[2]), sin(var_vals[3]))

G = [0, 0]
G[0] = groebner(poly_system[0], tuple(ccos[0:2]+ssin[0:4]+(sqrt3, ccos[2], ccos[3])), order='lex')
G[1] = groebner(poly_system[1], tuple(ccos[4:6]+ssin[2:6]+(sqrt3, ccos[2], ccos[3])), order='lex')

print("The polynomial system 1 implies (as can be proven via Groebner basis):")
pols_in_two_var = [g for g in G[1] if g.free_symbols <= {ccos[2], ccos[3]}]
assert(len(pols_in_two_var) == 1)
pol = pols_in_two_var[0]
print(pol.as_expr(), "= 0")
print("We can factor this polynomial and evaluate each factor in the case of the cubical configuration, which corresponds to c2 =", cos(var_vals[2]), "and c3 =", cos(var_vals[3]))
factors = factor(pol.as_expr()).as_ordered_factors()
for f in factors:
    eval_f = f.subs(ccos[2], cos(var_vals[2])).subs(ccos[3], cos(var_vals[3]))
    print(f, " -> ", eval_f.evalf())
    print()

print("Observe that all factors are nonzero apart from", factors[6], "which can be rewritten as", 1-(2*ccos[2]*ssin[3])**2, "=", (1-2*ccos[2]*ssin[3])*(1+2*ccos[2]*ssin[3]))
print("Only one of the two factors is zero in the cubical configuration, namely the first factor and thus we can add the polynomial", 2*ccos[2]*ssin[3]-1, "to polynomial system 0.")
print("By doing so, the polynomial system 0 (as can be proven via Groebner basis) implies that:", G[0][0].as_expr(), "is 0.")
print("From this point it is easy to conclude that the cubical configuration is rigid.")

print()
print("Let us remark that simply joining the two trigonometric/polynomial system does not yield the result effortlessly because the resulting algebraic system has a zero-set of dimension 1.")



# What follows is used to report the results in the paper.
# x, y = symbols('x y')
# display_pol = factor(pol.as_expr()).subs({ccos[2]: x, ccos[3]: y})
# print(display_pol, "= 0")