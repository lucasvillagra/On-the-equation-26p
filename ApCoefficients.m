/* Script to compute the forms with the right Fourier coefficients */

/* "d" is the coefficient of the equation, "new" the space of newforms, "j" is the number of the newform, "Primes" is the list of prime and "list" is the list of values of Chi evaluated in  "Primes", where the imput must be 4 if Chi equals i and -4 if equals -i.*/

ApCoefficients := function(d,new,j,Primes,list)

A:=[];
B:=[];
for p in Primes do
K:=Parent(Coefficient(Eigenform(new[j],12),11));
Comp:=AdjoinRoot(K,4);
Cyc<a>:=CyclotomicField(4);
bol1, map1:=IsSubfield(Cyc,Comp);
Root:=map1(a);
if KroneckerSymbol(-d,p) eq 1 then
A:=Append(A,Coefficient(Eigenform(new[j],90),p));
A:=Append(A,Coefficient(Eigenform(new[j],90),p));
else
A:=Append(A,Coefficient(Eigenform(new[j],90),p)^2-2*p*eps(p));
end if;
end for;
for n in [1..#A] do
if list[n] eq 4 then
B:=Append(B,A[n]*Root^-1);
end if;
if list[n] eq -4 then
B:=Append(B,A[n]/((-1)*Root));
end if;
if not IsDivisibleBy(list[n],4) then
B:=Append(B,A[n]/list[n]);
end if;
end for;
return(B);
end function;


/*======================================*/

CoefficientsCurve := function(E,Id)
B:=[];
for I in Id do
B:=Append(B,TraceOfFrobenius(E,I));
end for;
return(B);
end function;



