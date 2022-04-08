/* Case d=5 */

load "mazur.m";

print "Case d=5:";

/* Character chi */

Chi:= function(p)
Eps1:=Generators(DirichletGroup(15,CyclotomicField(4)));
Eps1:=Eps1[1]*Eps1[2];
Eps2:=Generators(DirichletGroup(20,CyclotomicField(4)));
Eps2:=Eps2[1]*Eps2[2];
if KroneckerSymbol(-5,p) eq -1 then
   f:=order(Eps1(p));
else
	f:=order(Eps2(p))*2;
end if;
return f;
end function;

/* First space */

print "Forms in Space 2^4*3^2*5^2:";
level1:="2^4*3^2*5^2";
G:=DirichletGroup(2^4*3^2*5^2,CyclotomicField(4));                           
eps:=(Elements(G)[50]);
M:=ModularSymbols(eps,2,1);
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
MZ:=DiscardPlace(5,eps,Chi,new,i,13,50);
print(MZ);
if MZ eq {@ 0 @} then
   BadForms1:=Append(BadForms1,i);
   end if;
end if;
end for;

print "Cannot discard the forms in the first space with parameter: ", BadForms1;

/* Second space */

print "Forms in Space 2^4*3^3*5^2:";
level2:="2^4*3^3*5^2";
eps2:=Extend(eps,2^4*3^3*5^2);
M2:=ModularSymbols(eps2,2,1);
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
MZ:=DiscardPlace(5,eps2,Chi,new2,i,13,50);
print(MZ);
if MZ eq {@ 0 @} then
   BadForms2:=Append(BadForms2,i);
   end if;
end if;
end for;

print "Cannot discard the forms in the second space with parameter: ", BadForms2;

BadForms:=[BadForms1,BadForms2];

/* Compute the curves attached to forms that cannot be discarded */

P<x> := PolynomialRing(Rationals());
K<a>:=NumberField(x^2+5);
O:=RingOfIntegers(K);
I:=ideal<O|30>;

/* Curves with good reduction out of primes dividing I */

S:=EllipticCurveWithGoodReductionSearch(I,700);

print "It was found", #S, "curves (over the", K, ") with good reduction out of primes dividing 2,3 or 5."; 

/* The primes to match forms and curves */

Id11:=ideal<O|11>;
Id13:=ideal<O|13>;
Id29_1:=ideal<O|3+2*a>;
Id29_2:=ideal<O|3-2*a>;
Id41_1:=ideal<O|6+a>;
Id41_2:=ideal<O|6-a>;
Id61_1:=ideal<O|4+3*a>;
Id61_2:=ideal<O|4-3*a>;
Id89_1:=ideal<O|-3+4*a>;
Id89_2:=ideal<O|-3-4*a>;

Cyc<a>:=CyclotomicField(4);

/* The character value at these ideals */

Char:=function(Id)
if Id eq Id11 then
   return -1;
   end if;
if Id eq Id13 then
   return -a;
   end if;
if Id eq Id29_1 then
   return -a;
   end if;
if Id eq Id29_2 then
   return a;
   end if;
if Id eq Id41_1 then
   return 1;
   end if;
if Id eq Id41_2 then
   return -1;
   end if;
if Id eq Id61_1 then
   return 1;
   end if;
if Id eq Id61_2 then
   return 1;
   end if;
if Id eq Id89_1 then
   return -a;
   end if;
if Id eq Id89_2 then
   return a;
   end if;
end function;   

/* Function to compute the coefficient of the base change */

Coef:=function(form,p)
if KroneckerSymbol(-5,p) eq 1 then
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

emb1:=[];
for g in new do
K:=Parent(Coefficient(Eigenform(g,12),11));
bol, map:=IsSubfield(Cyc,K);
emb1:=emb1 cat [map];
end for;

emb2:=[];
for g in new2 do
K:=Parent(Coefficient(Eigenform(g,12),11));
bol, map:=IsSubfield(Cyc,K);
emb2:=emb2 cat [map];
end for;

emb:=[emb1,emb2];

for E in S do
a11:=TraceOfFrobenius(E,Id11);
a13:=TraceOfFrobenius(E,Id13);
a29:=TraceOfFrobenius(E,Id29_1);
b29:=TraceOfFrobenius(E,Id29_2);
a41:=TraceOfFrobenius(E,Id41_1);
b41:=TraceOfFrobenius(E,Id41_2);
a61:=TraceOfFrobenius(E,Id61_1);
b61:=TraceOfFrobenius(E,Id61_2);
a89:=TraceOfFrobenius(E,Id89_1);
b89:=TraceOfFrobenius(E,Id89_2);

for j in [1..#BadForms] do
for i in BadForms[j] do
if [a11,a13,a29,b29,a41,b41,a61,b61,a89,b89] eq [Coef(New[j][i],11)*emb[j][i](Char(Id11))^3,Coef(New[j][i],13)*emb[j][i](Char(Id13))^3,Coef(New[j][i],29)*emb[j][i](Char(Id29_1))^3,Coef(New[j][i],29)*emb[j][i](Char(Id29_2))^3,Coef(New[j][i],41)*emb[j][i](Char(Id41_1))^3,Coef(New[j][i],41)*emb[j][i](Char(Id41_2))^3,Coef(New[j][i],61)*emb[j][i](Char(Id61_1))^3,Coef(New[j][i],61)*emb[j][i](Char(Id61_2))^3,Coef(New[j][i],89)*emb[j][i](Char(Id89_1))^3,Coef(New[j][i],89)*emb[j][i](Char(Id89_2))^3] then
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

Primes:=[11,13,29,41,61,89];
list:=[-1,-4,-4,4,1,-1,1,1,-4,4];


Id:=[Id11,Id13,Id29_1,Id29_2,Id41_1,Id41_2,Id61_1,Id61_2,Id89_1,Id89_2];


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
if ApCoefficients(5,New[j],i,Primes,list) eq C[j][k] then
A[j][k]:=Append(A[j][k],i);
end if;
end for;
print "List of forms of the Space", Level[j], "with Fourier coefficients a_p, with p in", Primes, "that match with the twist of each correspondent elliptic curve:", A[j][k];
end for;
end for;





