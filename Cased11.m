/* Case d=11 */

load "mazur.m";

print "Case d=11:";

/* Character chi */

Chi:= function(p)
if KroneckerSymbol(-11,p) eq -1 then
   f:=(3-KroneckerSymbol(p,33)) div 2;
else
	f:= 2;
end if;
return f;
end function;

/* eps (Nebentypus) */

eps:=Elements(DirichletGroup(1))[1];

/* First space */

print "Forms in Space 2^2*3^3*11^2:";
level1:="2^2*3^3*11^2";
M:=ModularSymbols(2^2*3^3*11^2,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);
CM:=FormsWithCM(new);
print "There are", #new, "forms,", #CM, "of them having complex multiplication.";
print "Forms with CM:", CM;

/* Mazur's trick for forms without CM */

print "Primes obtained via Mazur's trick for non-CM forms:";
BadForms1:=[];
for i in [1..#new] do         
if i notin CM then
MZ:=DiscardPlace(11,eps,Chi,new,i,13,50);
print(MZ);
if MZ eq {@ 0 @} then
   BadForms1:=Append(BadForms1,i);
   end if;
end if;
end for;

print "Cannot discard the forms in the first space with parameter: ", BadForms1;

/* Second space */

print "Forms in Space 2^2*3^2*11^2:";
level2:="2^2*3^2*11^2";
M2:=ModularSymbols(2^2*3^2*11^2,2,1);
S2:=NewSubspace(CuspidalSubspace(M2));
new2:=NewformDecomposition(S2);
CM2:=FormsWithCM(new2);
print "There are", #new2, "forms,", #CM2, "of them having complex multiplication.";
print "Forms with CM:", CM2;

Level:=[level1,level2];
New:=[new,new2];

/* Mazur's trick for forms without CM */

print "Primes obtained via Mazur's trick for non-CM forms:";
BadForms2:=[];
for i in [1..#new2] do 
if i notin CM2 then        
MZ:=DiscardPlace(11,eps,Chi,new2,i,5,40);
print(MZ);
if MZ eq {@ 0 @} then
   BadForms2:=Append(BadForms2,i);
   end if;
end if;
end for;

print "Cannot discard the forms in the second space with parameter: ", BadForms2;

BadForms:=[BadForms1, BadForms2];

/* Compute the curves attached to forms that cannot be discarded */

P<x> := PolynomialRing(Rationals());
K<a>:=NumberField(x^2+11);
O:=RingOfIntegers(K);
I:=ideal<O|66>;

/* Curves with good reduction out of primes dividing I */

S:=EllipticCurveWithGoodReductionSearch(I,700);

print "It was found", #S, "curves (over the", K, ") with good reduction out of primes dividing 2,3 or 11."; 

/* The primes to match forms and curves */

Id5_1:=ideal<O|(3+a)/2>;
Id5_2:=ideal<O|(3-a)/2>;
Id7:=ideal<O|7>;
Id13:=ideal<O|13>;
Id17:=ideal<O|17>;
Id19:=ideal<O|19>;
Id23_1:=ideal<O|(9+a)/2>;
Id23_2:=ideal<O|(-9+a)/2>;
Id31_1:=ideal<O|(5+3*a)/2>;
Id31_2:=ideal<O|(5-3*a)/2>;

Cyc<a>:=CyclotomicField(4);

/* The character value at these ideals */

Char:=function(Id)
if Id eq Id5_1 then
   return -1;
   end if;
if Id eq Id5_2 then
   return 1;
   end if;
if Id eq Id7 then
   return -1;
   end if;
if Id eq Id13 then
   return -1;
   end if;
if Id eq Id17 then
   return 1;
   end if;
if Id eq Id19 then
   return -1;
   end if;
if Id eq Id23_1 then
   return -1;
   end if;
if Id eq Id23_2 then
   return 1;
   end if;
if Id eq Id31_1 then
   return -1;
   end if;
if Id eq Id31_2 then
   return -1;
   end if;
end function;   


/* Function to compute the coefficient of the base change */

Coef:=function(form,p)
if KroneckerSymbol(-11,p) eq 1 then
return Coefficient(Eigenform(form,90),p);
else
return Coefficient(Eigenform(form,90),p)^2-2*p*eps(p);
end if;
end function;


/* The set Ans[j][i] will cointain all curves of S such that their ap's match the ap's of the newform i of the j-space */

Ans:=[];
for j in [1..#New] do
Ans[j]:=[];
for i in [1..#New[j]] do
Ans[j][i]:=[];
end for;
end for;

for E in S do
a5:=TraceOfFrobenius(E,Id5_1);
b5:=TraceOfFrobenius(E,Id5_2);
a7:=TraceOfFrobenius(E,Id7);
a13:=TraceOfFrobenius(E,Id13);
a17:=TraceOfFrobenius(E,Id17);
a19:=TraceOfFrobenius(E,Id19);
a23:=TraceOfFrobenius(E,Id23_1);
b23:=TraceOfFrobenius(E,Id23_2);
a31:=TraceOfFrobenius(E,Id31_1);
b31:=TraceOfFrobenius(E,Id31_2);

for j in [1..#BadForms] do
for i in BadForms[j] do
if [a5,b5,a7,a13,a17,a19,a23,b23,a31,b31] eq [Coef(New[j][i],5)*(Char(Id5_1)),Coef(New[j][i],5)*(Char(Id5_2)),Coef(New[j][i],7)*(Char(Id7)),Coef(New[j][i],13)*(Char(Id13)),Coef(New[j][i],17)*(Char(Id17)),Coef(New[j][i],19)*(Char(Id19)),Coef(New[j][i],23)*(Char(Id23_1)),Coef(New[j][i],23)*(Char(Id23_2)),Coef(New[j][i],31)*(Char(Id31_1)),Coef(New[j][i],31)*(Char(Id31_2))] then
Ans[j][i]:=Append(Ans[j][i],E);
end if;
end for;
end for;

end for;


for j in [1..#BadForms] do
for i in BadForms[j] do
print "List of curves corresponding to the newform", i,"of the Space", Level[j];
print(Ans[j][i]);
end for;
end for;

/* Forms attached to each Elliptic curve in Ans[j][i] */

load "ApCoefficients.m";

Primes:=[5,7,13,17,19,23,31];
list:=[-1,1,-1,-1,1,-1,-1,1,-1,-1];


Id:=[Id5_1,Id5_2,Id7,Id13,Id17,Id19,Id23_1,Id23_2,Id31_1,Id31_2];


C:=[];
for j in [1..#BadForms] do
C[j]:=[];
for i in BadForms[j] do
E:=Ans[j][i][1];
C[j]:=Append(C[j],CoefficientsCurve(E,Id));
end for;
end for;


A:=[];
for j in [1..#New] do
A[j]:=[];
for k in [1..#C[j]] do
A[j][k]:=[];
for i in [1..#New[j]] do
if ApCoefficients(11,New[j],i,Primes,list) eq C[j][k] then
A[j][k]:=Append(A[j][k],i);
end if;
end for;
print "List of forms of the Space", Level[j], "with Fourier coefficients a_p, with p in", Primes, "that match with the twist of each correspondent elliptic curve:", A[j][k];
end for;
end for;


