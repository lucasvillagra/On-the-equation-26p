// Muli-Frey function implemented just for the case d=4 and the particular newform that we must discard

MultiFrey:=function(new,p,Places,Order)
A:=[p];
Bad:=[];
if KroneckerSymbol(-4,p) eq 1 then
   Case:=1;
   else
   Case:=2;
end if;


G:=DirichletGroup(2^4*3);
eps:=Elements(G)[6];
eps:=Extend(eps,2^4*3^3);
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
New:=NewformDecomposition(S);
g:=New[3];  


K:=AbsoluteField(Parent(Coefficient(Eigenform(g,p+1),p)));
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

coefg:=Coefficient(Eigenform(g,p+1),p);



_<X> := PolynomialRing(Rationals());
L<r> := ext< Rationals() | X^2+4>;
R:=RingOfIntegers(L);
I:=ideal<R|p>;



for i in Places do
coeff:=Coefficient(Eigenform(new[i],p+1),p);




for a in [p..2*p] do
    for b in [p..2*p] do
            if not IsDivisibleBy(a*b,p^2) then
                if  not IsDivisibleBy(a^2+4*b^6,p) then 
                    coefE:=TraceOfFrobenius(EllipticCurve([6*b*r,0,-4*4*(a+b^3*r),0,0]),Factorization(I)[1][1]);
                    if Case eq 1 then
                          aux1:=Numerator(Norm(Norm(Norm(Norm(coefg-coefE*Root^pp))))); 
                          if  IsDivisibleBy(a^2+4*b^6,p) then
                          coefF:=TraceOfFrobenius(EllipticCurve([0,0,0,3*4*b^2,2*4*a]),p);
                          aux2:=Numerator(Norm(Norm(Norm(Norm(coefF-coeff)))));
                          else
                          aux2:=Numerator(Norm(coeff^2-(p+1)^2));
                          end if;
                    else
                          aux1:=Numerator(Norm(Norm(Norm(Norm(coefg^2-coefE*Root^pp-2*p*eps(p))))));
                          coefF:=TraceOfFrobenius(EllipticCurve([0,0,0,3*4*b^2,2*4*a]),p);
                          aux2:=Numerator(Norm(Norm(Norm(Norm(coefF-coeff)))));
                    end if;
                    
                    if Gcd(aux1,aux2) eq 0 then
                        Bad:=Append(Bad,i);
                    else
                        fact:=Factorisation(Gcd(aux1,aux2));
                           A:=A cat [g[1] : g in fact];
                    end if;
                    if Order ge 2 then
                        coefE:=TraceOfFrobenius(EllipticCurve([6*b*r,0,-4*4*(a+b^3*r),0,0]),Factorization(I)[1][1]);
                        if Case eq 1 then
                              aux1:=Numerator(Norm(Norm(Norm(Norm(coefg+coefE*Root^pp))))); 
                              if  IsDivisibleBy(a^2+4*b^6,p) then
                              coefF:=TraceOfFrobenius(EllipticCurve([0,0,0,3*4*b^2,2*4*a]),p);
                              aux2:=Numerator(Norm(Norm(Norm(Norm(coefF-coeff)))));
                              else
                                 aux2:=Numerator(Norm(coeff^2-(p+1)^2));
                              end if;                              
                        else
                              aux1:=Numerator(Norm(Norm(Norm(Norm(coefg^2+coefE*Root^pp-2*p*eps(p))))));
                              coefF:=TraceOfFrobenius(EllipticCurve([0,0,0,3*4*b^2,2*4*a]),p);
                              aux2:=Numerator(Norm(Norm(Norm(Norm(coefF-coeff)))));
                        end if;

                    if Gcd(aux1,aux2) eq 0 then
                        Bad:=Append(Bad,i);
                    else
                        fact:=Factorisation(Gcd(aux1,aux2));
                           A:=A cat [g[1] : g in fact];
                    end if;
                    end if;
                end if;
            end if;
    end for;
end for;
                    
end for;

if Case eq 1 then
A:=IndexedSet(A) join MazurTrickMultiplicative(New,p,eps,[3]);
end if;

return IndexedSet(A) , IndexedSet(Bad);
end function;


//======================================================================

DiscardPlaceF:=function(Chi,new,j,bd1,bd2)
st:=1;
for i in PrimesInInterval(bd1,bd2) do
if not (IsDivisibleBy(6,i)) then
   B,Bad:=MultiFrey(new,i,[j],Chi(i));
   if j notin Bad then
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
