/* Case d=15 */

load "mazur.m";

print "Case d=15:";

/* Character chi */

Chi:= function(p)
G:=DirichletGroup(2*3^5*5^2,CyclotomicField(4));
eps:=Elements(G)[4];
if KroneckerSymbol(-15,p) eq -1 then
   f:=order(eps(p)*KroneckerSymbol(p,5));
else
	f:=2*order(eps(p)*KroneckerSymbol(p,5));
end if;
return f;
end function;

/* First space */

print "Forms in Space 2*3^5*5^2:";
level1:="2*3^5*5^2";
G:=DirichletGroup(2*3^5*5^2,CyclotomicField(4));
eps:=Elements(G)[4];
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);
print "There are", #new, "forms";                       /* It is not possible to determine the CM forms with Magma in this case */

/* Mazur's trick for  all forms */

print "Primes obtained via Mazur's trick for all forms:";
BadForms1:=[];
for i in [1..#new-1] do         
MZ:=DiscardPlace(15,eps,Chi,new,i,5,40);
print(MZ);
if MZ eq {@ 0 @} then
   BadForms1:=Append(BadForms1,i);
   end if;
end for;

print(DiscardPlace(15,eps,Chi,new,#new,5,13));

if  DiscardPlace(15,eps,Chi,new,#new,5,13) eq {0} then
   BadForms1:=Append(BadForms1,#new);
end if;

print "Cannot discard the forms in the first space with parameter: ", BadForms1;

/* Second space */

print "Forms in Space 2^2*3^5*5^2:";
level2:="2^2*3^5*5^2";
G:=DirichletGroup(2^2*3^5*5^2,CyclotomicField(4));
eps2:=Elements(G)[7];
M2:=ModularSymbols(eps2,2,1);
S2:=NewSubspace(CuspidalSubspace(M2));
new2:=NewformDecomposition(S2);

print "There are", #new2, "forms";
print "The first form has complex multiplication:";
print "The form 1 has CM:", HasCM(new2[1]);
print "The forms 2,...,12 also have complex multiplication by -3 (see the paper)";

Level:=[level1,level2];
New:=[new,new2];

/* Mazur's trick for forms without CM */

print "Primes obtained via Mazur's trick for non-CM forms:";
BadForms2:=[];
for i in [13..23] do 
MZ:=DiscardPlace(15,eps,Chi,new2,i,5,40);
print(MZ);
if MZ eq {@ 0 @} then
   BadForms2:=Append(BadForms2,i);
   end if;
end for;

for i in [14..25] do 
MZ:=DiscardPlace(15,eps,Chi,new2,i,5,23);
print(MZ);
if MZ eq {@ 0 @} then
   BadForms2:=Append(BadForms2,i);
   end if;
end for;

for i in [16..#new] do 
MZ:=DiscardPlace(15,eps,Chi,new2,i,5,13);
print(MZ);
if MZ eq {@ 0 @} then
   BadForms2:=Append(BadForms2,i);
   end if;
end for;

print "Cannot discard the forms in the second space with parameter: ", BadForms2;

BadForms:=[BadForms1,BadForms2];

/* Compute the curves attached to forms that cannot be discarded */

P<x> := PolynomialRing(Rationals());
K<a>:=NumberField(x^2+15);
O:=RingOfIntegers(K);
I:=ideal<O|30>;

/* Curves with good reduction out of primes dividing I */

S:=EllipticCurveWithGoodReductionSearch(I,400);

print "It was found", #S, "curves (over the", K, ") with good reduction out of primes dividing 2,3 or 5."; 

/* The primes to match forms and curves */

Id7:=ideal<O|7>;
Id11:=ideal<O|11>;
Id13:=ideal<O|13>;
Id19_1:=ideal<O|2+a>;
Id19_2:=ideal<O|2-a>;
Id31_1:=ideal<O|4+a>;
Id31_2:=ideal<O|4-a>;
Id37:=ideal<O|37>;
Id41:=ideal<O|41>;

Cyc<a>:=CyclotomicField(4);

/* The character value at these ideals */

Char:=function(Id)
if Id eq Id7 then
   return -a;
   end if;
if Id eq Id11 then
   return -1;
   end if;
if Id eq Id13 then
   return a;
   end if;
if Id eq Id19_1 then
   return a;
   end if;
if Id eq Id19_2 then
   return a;
   end if;
if Id eq Id31_1 then
   return -1;
   end if;
if Id eq Id31_2 then
   return -1;
   end if;
if Id eq Id37 then
   return -a;
   end if;
if Id eq Id41 then
   return -1;
   end if;
end function;   


/* Function to compute the coefficient of the base change */

Coef:=function(form,p)
if KroneckerSymbol(-15,p) eq 1 then
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
Comp:=AdjoinRoot(K,4);
bol, map:=IsSubfield(Cyc,Comp);
emb1:=emb1 cat [map];
end for;

emb2:=[];
for g in new2 do
K:=Parent(Coefficient(Eigenform(g,12),11));
Comp:=AdjoinRoot(K,4);
bol, map:=IsSubfield(Cyc,Comp);
emb2:=emb2 cat [map];
end for;

emb:=[emb1,emb2];

for E in S do
a7:=TraceOfFrobenius(E,ideal<O|Id7>);
a11:=TraceOfFrobenius(E,ideal<O|Id11>);
a13:=TraceOfFrobenius(E,ideal<O|Id13>);
a19:=TraceOfFrobenius(E,ideal<O|Id19_1>);
b19:=TraceOfFrobenius(E,ideal<O|Id19_2>);
a31:=TraceOfFrobenius(E,ideal<O|Id31_1>);
b31:=TraceOfFrobenius(E,ideal<O|Id31_2>);
a37:=TraceOfFrobenius(E,ideal<O|Id37>);
a41:=TraceOfFrobenius(E,ideal<O|Id41>);

for j in [1..#BadForms] do
for i in BadForms[j] do
if [a7,a11,a13,a19,b19,a31,b31,a37,a41] eq [emb[j][i](Coef(New[j][i],7)*Char(Id7)^3),emb[j][i](Coef(New[j][i],11)*Char(Id11)^3),emb[j][i](Coef(New[j][i],13)*Char(Id13)^3),emb[j][i](Coef(New[j][i],19)*Char(Id19_1)^3),emb[j][i](Coef(New[j][i],19)*Char(Id19_2)^3),emb[j][i](Coef(New[j][i],31)*Char(Id31_1)^3),emb[j][i](Coef(New[j][i],31)*Char(Id31_2)^3),emb[j][i](Coef(New[j][i],37)*Char(Id37)^3),emb[j][i](Coef(New[j][i],41)*Char(Id41)^3)] then
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

Primes:=[7,11,13,19,31,37,41];
list:=[-4,-1,4,4,4,-1,-1,-4,-1];



Id:=[Id7,Id11,Id13,Id19_1,Id19_2,Id31_1,Id31_2,Id37,Id41];


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
if ApCoefficients(15,New[j],i,Primes,list) eq C[j][k] then
A[j][k]:=Append(A[j][k],i);
end if;
end for;
print "List of forms of the Space", Level[j], "with Fourier coefficients a_p, with p in", Primes, "that match with the twist of each correspondent elliptic curve:", A[j][k];
end for;
end for;


