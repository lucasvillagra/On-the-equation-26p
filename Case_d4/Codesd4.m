
/* Mazur's trick for case d=4, with the rational Frey curve F_(a,b,c) */


Mazur4:=function(new,p,Places)
A:=[p];
Bad:=[];

for i in Places do
   for a in [p..2*p] do
    for b in [p..2*p] do
        if not IsDivisibleBy(a*b,p^2) then
            if  not IsDivisibleBy(a^2+4*b^6,p) then
                apE:=TraceOfFrobenius(EllipticCurve([0,0,0,3*4*b^2,2*4*a]),p);
                aux:=Numerator(Norm(Coefficient(Eigenform(new[i],p+1),p)-apE));
            else
                aux:=Numerator(Norm(Coefficient(Eigenform(new[i],p+1),p)^2-(p+1)^2));
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
end for;

return IndexedSet(A), IndexedSet(Bad);
end function;


/* Mazur's trick for several primes */

DiscardPlace4:=function(new,j,bd1,bd2)
st:=1;
for i in PrimesInInterval(bd1,bd2) do
if not (IsDivisibleBy(6,i)) then
   B,Bad:=Mazur4(new,i,[j]);
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
