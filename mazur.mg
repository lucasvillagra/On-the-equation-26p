/* Scripts needed to discard forms. In all cases, new is a list of newforms obtained by NewformDecomposition. */

/* Script to compute all forms with CM in the list */

FormsWithCM:=function(new)
A:=[];
for i in [1 .. #new] do
if HasCM(new[i]) then
A:=Append(A,i);
end if;
end for;
return A;
end function;


/* Script to given d and p, compute the possitble values (a,b) of the equation modulo p, and the ap of each such curve. */


ApCand:=function(d,p)
list:=[];
_<X> := PolynomialRing(Rationals());
K<r> := ext< Rationals() | X^2+d>;
R:=RingOfIntegers(K);
I:=ideal<R|p>;
for a in [p..2*p] do
    for b in [p..2*p] do
        if not IsDivisibleBy(a*b,p^2) then
            if  not IsDivisibleBy(a^2+d*b^6,p) then
                list:=Append(list,TraceOfFrobenius(EllipticCurve([6*b*r,0,-4*d*(a+b^3*r),0,0]),Factorization(I)[1][1]));
            end if;
        end if;
    end for;
end for;
list:=Set(list);
return(list);
end function;


/* This are internal routines used in the previous script. The first one is needed to compute values of the character chi on the coefficient field of the newform. */

AdjoinRoot := function(K, p: DoLinearExtension := false)
  f := CyclotomicPolynomial(p);
  f := Polynomial(K, f);
  f := Factorisation(f)[1][1];
  return AbsoluteField(ext<K|f : DoLinearExtension := DoLinearExtension>);
end function;


/* This script corresponds to the first two cases of Mazur's trick (i.e. do not consider the lowering the level situation). Input: d, newformsdecomposition, prime p, Nebentypus, the places of new to discard, the elliptic curve ap's, and the order of chi at the prime p. */

MazurTrick:=function(d,new,p,eps,Places,ApCand,Order)
A:=[p];
Bad:=[];
if KroneckerSymbol(-d,p) eq 1 then
   Case:=1;
   else
   Case:=2;
end if;
for i in Places do
K:=AbsoluteField(Parent(Coefficient(Eigenform(new[i],p+1),p)));
if K eq Rationals() then
   K:=RationalsAsNumberField();
   end if;
Comp:=AdjoinRoot(K,Order);
Cyc<a>:=CyclotomicField(Order);
bol1, map1:=IsSubfield(Cyc,Comp);
Root:=map1(a);
pp:=1;
while Root^(2*pp) ne eps(p^Case) do
   pp:= pp+2;
   end while;
for j in ApCand do                                                            
   if Case eq 1 then
      aux:=Numerator(Norm(Norm(Norm(Norm(Coefficient(Eigenform(new[i],p+1),p)-j*Root^pp)))));
      else
      aux:=Numerator(Norm(Norm(Norm(Norm(Coefficient(Eigenform(new[i],p+1),p)^2-j*Root^pp-2*p*eps(p))))));
   end if;
   if aux eq 0 then
      Bad:=Append(Bad,i);
      else
      fact:=Factorisation(aux);
      A:=A cat [g[1] : g in fact];
      end if;
if Order ge 2 then
   if Case eq 1 then
      aux:=Numerator(Norm(Norm(Norm(Norm(Coefficient(Eigenform(new[i],p+1),p)+j*Root^pp)))));
      else
      aux:=Numerator(Norm(Norm(Norm(Norm(Coefficient(Eigenform(new[i],p+1),p)^2+j*Root^pp-2*p*eps(p))))));
   end if;
   if aux eq 0 then
      Bad:=Append(Bad,i);
      else
      fact:=Factorisation(aux);
      A:=A cat [g[1] : g in fact];
      end if;
end if;
end for;
end for;
return IndexedSet(A), IndexedSet(Bad);
end function;


/* This is Mazur's trick in the lowering the level case. */

MazurTrickMultiplicative:=function(new,p,eps,Places)
A:=[];
for i in Places do
aux:=Numerator(Norm(Norm(Coefficient(Eigenform(new[i],p+1),p)^2-1/eps(p)*(p+1)^2)));
fact:=Factorisation(aux);
A:=[g[1] : g in fact] cat A;
end for;
return IndexedSet(A);
end function;

/* This is the main script used to discard the form in the place j of the newformdecomposition. 
The bounds bd1 and bd2 give the range of primes used in Mazur's trick. */

DiscardPlace:=function(d,eps,Chi,new,j,bd1,bd2)
st:=1;
for i in PrimesInInterval(bd1,bd2) do
if not (IsDivisibleBy(6*d,i)) then
   B,Bad:=MazurTrick(d,new,i,eps,[j],ApCand(d,i),Chi(i));
   if j notin Bad then
   if KroneckerSymbol(-d,i) eq 1 then
        B:= B join MazurTrickMultiplicative(new,i,eps,[j]);
   end if;
   if st eq 1 then
      C:=B;
      st:=0;
      else 
      C:=C meet B;
      end if;
   end if;
end if;
end for;
if st eq 0 then
return C;
else
print "Cannot discard the form with parameter: ", j;
return({0});
end if;
end function;

/* This script will be used to define the character chi in each case */

order:=function(a)
if a eq 1 then
   f:=1;
else
   if a^2 eq 1 then
      f:=2;
      else
      f:=4;
      end if;
   end if;
return f;
end function;
