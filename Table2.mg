/* Script to compute Table 2.1 */

/* Conductor at 2 (for d square-free) */

cond2:= function(d)
    if IsDivisibleBy(d,2) then
        C:={2,3,4,7};
    else
        C:={0,1,2,3,4,5,6};
    end if;
return C;
end function;

/* Conductor at 3 (for d square-free) */

cond3:= function(d)
    if IsDivisibleBy(d,3) then
        C:={5};
    else
        C:={2,3};
    end if;
return C;
end function;

/* Function to compute Table 2.1 (for d \neq 19) */

ConductorCurve:=function(d)
D:= CremonaDatabase();
minC, maxC := ConductorRange(D);
dt:= d div (2^Valuation(d,2)) div (3^Valuation(d,3));
for a in cond2(d)  do for b in cond3(d) do 
  N:=2^a*3^b*dt^2;
  F:=EllipticCurves(D,N);
  for f in F do 
    c:=cInvariants(f);
    if (IsDivisibleBy(c[1],d)) and (IsDivisibleBy(c[2],d)) and (c[1] lt 0) and (c[2] lt 0) and (IsDivisibleBy(c[1],9)) and (IsDivisibleBy(Round(c[2]),27)) and (Discriminant(f) lt 0) then 
  T:=Round(-c[1]/d);
  if IsSquare(T) then 
    if (IsDivisibleBy(Round(c[1]),16)) and (IsDivisibleBy(Round(c[2]),64)) 
      then y:=Round(Sqrt(Round(-c[1]/(144*d)))); x:=Round(-c[2]/(1728*d)); 
         if Gcd(x,d*y) eq 1 then 
            print f ," of conductor", N, ". In this case [a,b]=", [x,y], "and p is a divisor of", Gcd(Valuation(x^2+d*y^6,2),Valuation(x^2+d*y^6,3));
         end if;
    end if;
    if ((Round(c[1]) mod 2) eq 1) and ((Round(c[2]) mod 2) eq 1) 
      then y:=Round(Sqrt(Round(-c[1]/(9*d)))); x:=Round(-c[2]/(27*d)); 
           if Gcd(x,d*y) eq 1 then 
                print f, ", of conductor", N, ". In this case [a,b]=", [x,y], "and p is a divisor of", Gcd(Valuation(x^2+d*y^6,2),Valuation(x^2+d*y^6,3));
           end if;
    end if;
end if;
end if; 
end for; end for; end for;
return {};
end function;
