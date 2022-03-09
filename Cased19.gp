\\ Case d=19, for 16 (see "Case19.mg")

print("Case d=19, discarding form 16:")

P=read("pol5.txt");
coef5=read("coef5.txt");
coef11=read("coef11.txt");
P=subst(P,x,a);

n5=norm(coef5*Mod(1,P))*norm(coef5^2*Mod(1,P)+9);
n11=norm(coef11*Mod(1,P))*norm(coef11^2*Mod(1,P)+9)*norm(coef11^2*Mod(1,P)+36);

print("Primes obtained via Mazur's trick for form 16:")
factor(gcd(n5,n11))
