/* Case d=13 */

print "Case d=13:";

load "mazur.m";

/* Character chi */

/* This function just compure the order of Chi. We use it to run Mazur's trick as Section 1.1. It will be an imput in the function "DiscardPlace" */

Chi:=function(p)
if KroneckerSymbol(-13,p) eq -1 then
    f:=1;
else
    f:=2*order(KroneckerSymbol(-1,p)*KroneckerSymbol(p,3)*KroneckerSymbol(p,13));
end if;
return f;
end function;

/* First space */

print "Forms in Space 2^4*3*13^2:";
G:=DirichletGroup(2^4*3*13^2);
eps:=Elements(G)[14];
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);
CM:=FormsWithCM(new);
print "There are", #new, "forms,", #CM, "of them having complex multiplication.";
print "Forms with CM:", CM;

/* Mazur's trick for forms without CM */

print "Primes obtained via Mazur's trick for non-CM forms:";
for i in [1..27] do         
if i notin CM then
print(DiscardPlace(13,eps,Chi,new,i,5,40));
end if;
end for; 

for i in [28..#new] do         
if i notin CM then
print(DiscardPlace(13,eps,Chi,new,i,5,10));
end if;
end for;

/* Second space */

print "Forms in Space 2^4*3^3*13^2:";
G:=DirichletGroup(2^4*3^3*13^2);
eps:=Elements(G)[14];
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);
CM:=FormsWithCM(new);
print "There are", #new, "forms,", #CM, "of them having complex multiplication.";
print "Forms with CM:", CM;

/* Mazur's trick for forms without CM */

print "Primes obtained via Mazur's trick for non-CM forms:";
for i in [1..57] do         
if i notin CM then
print(DiscardPlace(13,eps,Chi,new,i,5,40));
end if;
end for;

for i in [58..#new] do         
if i notin CM then
print(DiscardPlace(13,eps,Chi,new,i,5,9));
end if;
end for; 
