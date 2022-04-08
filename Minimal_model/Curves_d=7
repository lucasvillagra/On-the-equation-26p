# Case d=7 (Sage)

K.<a> = NumberField(x^2+7)


# Space 2*3*7^2

# NEWFORM 1


E=EllipticCurve([1,2,1/2*(a - 1),(950*a - 44),(4593*a + 100179)])
E.global_minimal_model()

# Output: Elliptic Curve defined by y^2 + x*y + (1/2*a+1/2)*y = x^3 + (-1)*x^2 + (950*a-46)*x + (7285/2*a+200449/2) over Number Field in a with defining polynomial x^2 + 7

# ==================================================================================================

# Space 2*3^3*7^2

# NEWFORM1

E=EllipticCurve([1,2,1/2*(a - 1),1/2*(1039*a - 823),-5977*a - 7131])
E=E.global_minimal_model(semi_global=True)

# Output: Elliptic Curve defined by y^2 + x*y + (1/2*a+1/2)*y = x^3 + (-1)*x^2 + (1039/2*a-827/2)*x + (-6497*a-6718) over Number Field in a with defining polynomial x^2 + 7

F=EllipticCurve([1,-1,1/2*(-1 + a),1/2*(115*a-91),1/2*(443*a+529)])

F.is_quadratic_twist(E)
# Output: -1/3

# ==================================================================================================

# NEWFORM 2


E=EllipticCurve([1,-1,1/2*(-a + 1),(16*a -46),1/2*(121*a - 157)])
E.global_minimal_model()

# Output: Elliptic Curve defined by y^2 + x*y + (1/2*a-1/2)*y = x^3 + (-1)*x^2 + (31/2*a-91/2)*x + (121/2*a-157/2) over Number Field in a with defining polynomial x^2 + 7

# ==================================================================================================

# NEWFORM3


E=EllipticCurve([-1,-1,1/2*(a - 1),(-26*a -46),1/2*(-201*a - 59)])
E=E.global_minimal_model()

# Output: Elliptic Curve defined by y^2 + x*y + (1/2*a-1/2)*y = x^3 + (-1)*x^2 + (-53/2*a-91/2)*x + (-201/2*a-59/2) over Number Field in a with defining polynomial x^2 + 7


F=EllipticCurve([1,-1,1/2*(-1+a),1/2*(-91-53*a),-1/2*(201*a+59)])

E.is_isogenous(F)
# Output: True
E.is_quadratic_twist(F)
# Output: 1
