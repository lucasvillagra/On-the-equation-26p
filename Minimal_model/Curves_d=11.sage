# Case d=11 (Sage)

K.<a> = NumberField(x^2+11)

# Space 2^2*3^2*11^2

# NEWFORM 15

E=EllipticCurve([0,-2*a-4,0,-948*a - 816,21312*a - 35028])
E=E.global_minimal_model()
E
# Output: Elliptic Curve defined by y^2 = x^3 + (-1/2*a+1/2)*x^2 + (-1907/2*a-1615/2)*x + (19479*a-31012) over Number Field in a with defining polynomial x^2 + 11

F=EllipticCurve([0,(1/2*a+1/2),0,(1907/2*a-1615/2),(-19479*a-31012)])

E.is_isogenous(F.quadratic_twist(-3))
# Output: True
