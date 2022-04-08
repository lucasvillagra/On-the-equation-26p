# Case d=5 (Sage)

K.<a> = NumberField(x^2+5)

# Space 2^4*3^2*5^2


# NEWFORM8

R.<x> = QQ[]
K.<a> = NumberField(x^2+5)
E=EllipticCurve([0,0,0,(855360*a + 388800),(212284800*a + 919123200)])
E.global_minimal_model(semi_global=True)

# Output: Elliptic Curve defined by y^2 + y = x^3 + (a-1)*x^2 + (a-43)*x + (-21*a+113) over Number Field in a with defining polynomial x^2 + 5

# ============================================================================================

# NEWFORM11


E=EllipticCurve([0,0,0,(-349920*a + 2843100),(-716461200*a - 724334400)])
E.global_minimal_model(semi_global=True)

# Output: Elliptic Curve defined by y^2 + (a+1)*x*y + (a+1)*y = x^3 + (-a)*x^2 + (-4*a+29)*x + (-25*a-32) over Number Field in a with defining polynomial x^2 + 5


# ============================================================================================

# NEWFORM12

E=EllipticCurve([0,0,0,(-172800*a - 418500),(67662000*a+64368000)])
E.global_minimal_model(semi_global=True)

# Output: Elliptic Curve defined by y^2 + (a+1)*x*y = x^3 + (-1)*x^2 + (-6*a-12)*x + (-10*a+2) over Number Field in a with defining polynomial x^2 + 5

# ---------------------------------------------------------------------

E=EllipticCurve([0,0,0,(1555200*a - 486000),(-653184000*a - 1737936000)])
E.global_minimal_model(semi_global=True)

# Output: Elliptic Curve defined by y^2 + (a+1)*x*y = x^3 + (-1)*x^2 + (-a-2)*x + 4 over Number Field in a with defining polynomial x^2 + 5

# ============================================================================================

# NEWFORM13


E=EllipticCurve([0,0,0,(100932480*a + 122860800),(55781913600*a + 1362341203200)])   
E.global_minimal_model(semi_global=True)

# Output: Elliptic Curve defined by y^2 = x^3 + (-a+1)*x^2 + (-1054*a-668)*x + (-16514*a+14618) over Number Field in a with defining polynomial x^2 + 5


# ============================================================================================

# Space 2^4*3^3*5^2

# NEWFORM 7


E=EllipticCurve([0,0,0,(-2138400*a + 3402000),(-408240000*a - 5534568000)])   
E.global_minimal_model(semi_global=True)

# Output: Elliptic Curve defined by y^2 = x^3 + (24*a+60)*x + (-56*a+256) over Number Field in a with defining polynomial x^2 + 5


# ============================================================================================

# NEWFORM8

E=EllipticCurve([0,0,0,(-58320*a - 72900),(-10206000*a + 583200)])   
E.global_minimal_model()

# Output: Elliptic Curve defined by y^2 + (a+1)*x*y + (a+1)*y = x^3 + (a+1)*x^2 + (-46*a-58)*x + (-248*a+126) over Number Field in a with defining polynomial x^2 + 5

