/* Case d=7 */

load "mazur.m";

print "Case d=7:";

/* Character chi */

Chi:= function(p)
if KroneckerSymbol(-7,p) eq -1 then
   f:=1;
else
	f:=KroneckerSymbol(p,3)*KroneckerSymbol(p,7);
    if f eq 1 then
        f:=2;
    else 
        f:=4;
    end if;
end if;
return f;
end function;

/* Fist space */

print "Forms in Space 2*3*7^2:";
level1:="2*3*7^2";
G:=DirichletGroup(2^1*3^1*7^2);
eps:=Generators(G)[1]*Generators(G)[2];
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
MZ:=DiscardPlace(7,eps,Chi,new,i,13,23);
print(MZ);
if MZ eq { 0 } then
   BadForms1:=Append(BadForms1,i);
   end if;
end if;
end for;

print "Cannot discard the forms in the first space with parameter: ", BadForms1;

/* Second space */

print "Forms in Space 2^2*3*7^2:";
level2:="2^2*3*7^2";
G:=DirichletGroup(2^2*3*7^2);
eps2:=Generators(G)[2]*Generators(G)[3];
M:=ModularSymbols(eps2,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new2:=NewformDecomposition(S);
CM2:=FormsWithCM(new2);

print "There are", #new2, "forms,", #CM2, "of them having complex multiplication.";
print "Forms with CM:", CM2;


/* Mazur's trick for forms without CM */

print "Primes obtained via Mazur's trick for non-CM forms:";
BadForms2:=[];
for i in [1..#new2] do 
if i notin CM2 then        
MZ:=DiscardPlace(7,eps,Chi,new2,i,2,20);
print(MZ);
if MZ eq {@ 0 @} then
   BadForms2:=Append(BadForms2,i);
   end if;
end if;
end for;

print "Cannot discard the forms in the second space with parameter: ", BadForms2;

/* Third space */

print "Forms in Space 2^2*3^3*7^2:";
level3:="2^2*3^3*7^2";
G:=DirichletGroup(2^2*3^3*7^2);
eps3:=Generators(G)[2]*Generators(G)[3];
M:=ModularSymbols(eps3,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new3:=NewformDecomposition(S);
CM3:=FormsWithCM(new3);

print "There are", #new3, "forms,", #CM3, "of them having complex multiplication.";
print "Forms with CM:", CM3;

/* Mazur's trick for forms without CM */

print "Primes obtained via Mazur's trick for non-CM forms:";
BadForms3:=[];
for i in [1..#new3] do 
if i notin CM3 then        
MZ:=DiscardPlace(7,eps,Chi,new3,i,2,20);
print(MZ);
if MZ eq {@ 0 @} then
   BadForms3:=Append(BadForms3,i);
   end if;
end if;
end for;

print "Cannot discard the forms in the third space with parameter: ", BadForms3;

/* Fourth space */

print "Forms in Space 2*3^3*7^2:";
level4:="2*3^3*7^2";
G:=DirichletGroup(2*3^3*7^2);
eps4:=Generators(G)[1]*Generators(G)[2];
M:=ModularSymbols(eps4,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new4:=NewformDecomposition(S);
CM4:=FormsWithCM(new4);

print "There are", #new4, "forms,", #CM4, "of them having complex multiplication.";
print "Forms with CM:", CM4;

Level:=[level1,level2,level3,level4];
New:=[new,new2,new3,new4];

/* Mazur's trick for forms without CM */

print "Primes obtained via Mazur's trick for non-CM forms:";
BadForms4:=[];
for i in [1..#new4] do 
if i notin CM4 then        
MZ:=DiscardPlace(7,eps,Chi,new4,i,2,13);
print(MZ);
if MZ eq {@ 0 @} then
   BadForms4:=Append(BadForms4,i);
   end if;
end if;
end for;

print "Cannot discard the forms in the fourth space with parameter: ", BadForms4;

BadForms:=[BadForms1,BadForms2, BadForms3, BadForms4];

/* Compute the curves attached to forms that cannot be discarded */

P<x> := PolynomialRing(Rationals());
K<a>:=NumberField(x^2+7);
O:=RingOfIntegers(K);
I:=ideal<O|42>;

/* Curves with good reduction out of primes dividing I */

S:=EllipticCurveWithGoodReductionSearch(I,700);

print "It was found", #S, "curves (over the", K, ") with good reduction out of primes dividing 2,3 or 7."; 

/* The primes to match forms and curves */

Id5:=ideal<O|5>;
Id11_1:=ideal<O|2+a>;
Id11_2:=ideal<O|2-a>;
Id13:=ideal<O|13>;
Id17:=ideal<O|17>;
Id19:=ideal<O|19>;
Id23_1:=ideal<O|4+a>;
Id23_2:=ideal<O|4-a>;
Id29_1:=ideal<O|-1+2*a>;
Id29_2:=ideal<O|-1-2*a>;
Id37_1:=ideal<O|3+2*a>;
Id37_2:=ideal<O|3-2*a>;
Id53_1:=ideal<O|-5+2*a>;
Id53_2:=ideal<O|-5-2*a>;


Cyc<a>:=CyclotomicField(4);

/* The character value at these ideals */

Char:=function(Id)
if Id eq Id5 then
   return 1;
   end if;
if Id eq Id11_1 then
   return -a;
   end if;
if Id eq Id11_2 then
   return a;
   end if;
if Id eq Id13 then
   return 1;
   end if;
if Id eq Id17 then
   return 1;
   end if;
if Id eq Id19 then
   return 1;
   end if;
if Id eq Id23_1 then
   return a;
   end if;
if Id eq Id23_2 then
   return -a;
   end if;
if Id eq Id29_1 then
   return a;
   end if;
if Id eq Id29_2 then
   return -a;
   end if;
if Id eq Id37_1 then
   return -1;
   end if;
if Id eq Id37_2 then
   return -1;
   end if;
if Id eq Id53_1 then
   return -a;
   end if;
if Id eq Id53_2 then
   return a;
   end if;

end function;   

/* Function to compute the coefficient of the base change */

Coef:=function(form,p)
if KroneckerSymbol(-7,p) eq 1 then
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

emb3:=[];
for g in new3 do
K:=Parent(Coefficient(Eigenform(g,12),11));
Comp:=AdjoinRoot(K,4);
bol, map:=IsSubfield(Cyc,Comp);
emb3:=emb3 cat [map];
end for;

emb4:=[];
for g in new4 do
K:=Parent(Coefficient(Eigenform(g,12),11));
Comp:=AdjoinRoot(K,4);
bol, map:=IsSubfield(Cyc,Comp);
emb4:=emb4 cat [map];
end for;

emb:=[emb1, emb2, emb3, emb4];


for E in S do
a5:=TraceOfFrobenius(E,Id5);
a11:=TraceOfFrobenius(E,Id11_1);
b11:=TraceOfFrobenius(E,Id11_2);
a13:=TraceOfFrobenius(E,Id13);
a17:=TraceOfFrobenius(E,Id17);
a19:=TraceOfFrobenius(E,Id19);
a23:=TraceOfFrobenius(E,Id23_1);
b23:=TraceOfFrobenius(E,Id23_2);
a29:=TraceOfFrobenius(E,Id29_1);
b29:=TraceOfFrobenius(E,Id29_2);
a37:=TraceOfFrobenius(E,Id37_1);
b37:=TraceOfFrobenius(E,Id37_2);
a53:=TraceOfFrobenius(E,Id53_1);
b53:=TraceOfFrobenius(E,Id53_2);

for j in [1..#BadForms] do
for i in BadForms[j] do
if [a5,a11,b11,a13,a17,a19,a23,b23,a29,b29,a37,b37,a53,b53] eq [emb[j][i](Coef(New[j][i],5)*Char(Id5)^3),emb[j][i](Coef(New[j][i],11)*Char(Id11_1)^3),emb[j][i](Coef(New[j][i],11)*Char(Id11_2)^3),emb[j][i](Coef(New[j][i],13)*Char(Id13)^3),emb[j][i](Coef(New[j][i],17)*Char(Id17)^3),emb[j][i](Coef(New[j][i],19)*Char(Id19)^3),emb[j][i](Coef(New[j][i],23)*Char(Id23_1)^3),emb[j][i](Coef(New[j][i],23)*Char(Id23_2)^3),emb[j][i](Coef(New[j][i],29)*Char(Id29_1)^3),emb[j][i](Coef(New[j][i],29)*Char(Id29_2)^3),emb[j][i](Coef(New[j][i],37)*Char(Id37_1)^3),emb[j][i](Coef(New[j][i],37)*Char(Id37_2)^3),emb[j][i](Coef(New[j][i],53)*Char(Id53_1)^3),emb[j][i](Coef(New[j][i],53)*Char(Id53_2)^3)] then
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

Primes:=[5,11,13,17,19,23,29,37,53];
list:=[1,-4,4,1,1,1,4,-4,4,-4,-1,-1,-4,4];


Id:=[Id5,Id11_1,Id11_2,Id13,Id17,Id19,Id23_1,Id23_2,Id29_1,Id29_2,Id37_1,Id37_2,Id53_1,Id53_2];


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
if ApCoefficients(7,New[j],i,Primes,list) eq C[j][k] then
A[j][k]:=Append(A[j][k],i);
end if;
end for;
print "List of forms of the Space", Level[j], "with Fourier coefficients a_p, with p in", Primes, "that match with the twist of each correspondent elliptic curve:", A[j][k];
end for;
end for;


