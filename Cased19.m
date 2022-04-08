/* Case d=19 */

load "mazur.mg";

print "Case d=19:";

/* Character chi */

Chi:= function(p)
if KroneckerSymbol(-19,p) eq -1 then
   f:=1;
else
	f:=2*order(KroneckerSymbol(p,3)*KroneckerSymbol(p,19));
end if;
return f;
end function;

/* First space */

print "Forms in Space 2^2*3*19^2:";
G:=DirichletGroup(2^2*3*19^2);
eps:=Elements(G)[7];
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);
CM:=FormsWithCM(new);
print "There are", #new, "forms,", #CM, "of them having complex multiplication.";
print "Forms with CM:", CM;

/* Mazur's trick for forms without CM */

print "Primes obtained via Mazur's trick for non-CM forms:";
for i in [1..#new] do
if i notin CM then
print(DiscardPlace(19,eps,Chi,new,i,2,30));
end if;
end for;

/* Second space */

print "Forms in Space 2^2*3^3*19^2:";
G:=DirichletGroup(2^2*3^3*19^2);
eps:=Elements(G)[7];
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);

/* Three forms with CM */
print "The forms 1, 2 and 7 have CM:";
for i in [1,2,7] do
print "HasCM(new[",i,"]));";
print(HasCM(new[i]));
end for;

/* Mazur's trick for all forms but 1,2,7,16 */

print "Primes obtained via Mazur's trick for all forms but 1,2,7,16:";
for i in [3..5] do
print(DiscardPlace(19,eps,Chi,new,i,2,50));
end for;

print(DiscardPlace(19,eps,Chi,new,6,2,70));

for i in [8..15] do
print(DiscardPlace(19,eps,Chi,new,i,2,50));
end for;

for i in [17..#new] do
print(DiscardPlace(19,eps,Chi,new,i,2,13));
end for;

print "The form 16 is discarded using PARI/GP (see 'Case19.gp').";

/* Discarding new[16] */

/* Computing coefficients a5 and a11 to discard form 16 */

f:=PowerSeries(new[16],40);
Write("coef5.txt",Coefficient(f,5));
Write("coef11.txt",Coefficient(f,11));
Write("pol5.txt",Parent(Coefficient(f,5)));

/* Computing possible values of a5 and a11 for the Frey curve */

A5:=ApCand(19,5);

A11:=ApCand(19,11);

/* See "Case19.gp" */
