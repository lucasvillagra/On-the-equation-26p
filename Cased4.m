/* Case d=4 */

load "mazur.m";

load "Multi-Frey.m";

print "Case d=4:";

/* Character chi */

/* This function just compure the order of Chi. We use it to run Mazur's trick as Section 1.1. It will be an imput in the function "DiscardPlace" */

Chi:= function(p)
G:=DirichletGroup(2^4*3);
eps:=Elements(G)[6];
if KroneckerSymbol(-1,p) eq -1 then
   f:=order(eps(p));
else
	f:=order(eps(p))*2;
end if;
return f;
end function;

/* First space */

print "Forms in Space 2^4*3^2*5^2:";
level1:="2^4*3";
G:=DirichletGroup(2^4*3);
eps:=Elements(G)[6];
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
MZ:=DiscardPlace(4,eps,Chi,new,i,5,30);
print(MZ);
if MZ eq {@ 0 @} then
   BadForms1:=Append(BadForms1,i);
   end if;
end if;
end for;

print "Cannot discard the forms in the first space with parameter: ", BadForms1;

/* Second space */

print "Forms in Space 2^4*3^3*5^2:";
level2:="2^4*3^3";
eps2:=Extend(eps,2^4*3^3);
M2:=ModularSymbols(eps2,2,1);
S2:=NewSubspace(CuspidalSubspace(M2));
new2:=NewformDecomposition(S2);
CM2:=FormsWithCM(new2);
print "There are", #new2, "forms,", #CM2, "of them having complex multiplication.";
print "Forms with CM:", CM2;


/* Mazur's trick for forms without CM */

print "Primes obtained via Mazur's trick for non-CM forms:";
BadForms2:=[];
for i in [1..#new2] do 
if i notin CM2 then        
MZ:=DiscardPlace(4,eps2,Chi,new2,i,3,50);
print(MZ);
if MZ eq {@ 0 @} then
   BadForms2:=Append(BadForms2,i);
   end if;
end if;
end for;

print "Cannot discard the forms in the second space with parameter: ", BadForms2;

/* Multi-Frey  */

print "Primes obtained via Multi-Frey to discard the third newform of the space of level 2^3*3^3:";

/* While comparing with the space of level 2^6*3, associated to the rational Frey elliptic curve */

G:=DirichletGroup(2^6*3);
eps:=Elements(G)[1];
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);


Bad:=[];
for i in [1..#new] do         
MZ:=DiscardPlaceF(Chi,new,i,5,11);
print(MZ);
if MZ eq {@ 0 @} then
   Bad:=Append(Bad,i);
end if;
end for;


/* While comparing with the space of level 2^6*3^3, associated to the rational Frey elliptic curve */

eps:=Extend(eps,2^6*3^3);
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);


for i in [1..#new] do         
MZ:=DiscardPlaceF(Chi,new,i,5,11);
print(MZ);
if MZ eq {@ 0 @} then
   Bad:=Append(Bad,i);
end if;
end for;

assert Bad eq [];  // This implies that the third newform of the space of level 2^3*3^3 is not more a "bad" newform
