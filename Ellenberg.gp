\\======================================================================
\\ Ellenberg. q = conductor character attached to K/Q, following Bennett-Ellenberg


\\ This function corresponds to the first bound in Proposition 9 (N>400). The variable s = sigma in their notation.
\\ Take s=q^2/2/Pi;

bound1(m,N,q,s)=2*(4*Pi*zeta(3/2)^2+1)*(400/399)^3*exp(2*Pi)*q^2*m^(3/2)*sigma(N,0)/N^(1/2+2*Pi*s/q^2);

E1(m,N,q,s)=16/3*Pi^3*m^(3/2)*s*log(N)/exp(N/2/Pi/m/s/log(N));

E3(m,N,q,s)=8/3*zeta(3/2)^2*Pi^3*s*m^(3/2)*sigma(N,0)*log(N)/sqrt(N)/exp(N/2/Pi/m/s/log(N));


\\ The next one comes from Proposition 10. Assume q is even, N>400.

E2(m,N,q,s)=
{local(o,x,Z);
if(s,,s=q^2/2/Pi);
Z=lfuncreate(1);
x =s*N*log(N);
o=exp(-2*Pi/x);
if(q%2==0,64*q*eulerphi(q)*Pi^5*m^2*(zeta(2)/6/N^2+1/Pi*(zeta(3)*log(exp(1)*N/2)-lfun(Z,3,1))/N^3)+32*Pi^5*zeta(7/2)^2*m^(5/2)*sigma(N,0)/N^(7/2)*((N^2/4/Pi^2/m+1)/(1-o)+1/(1-o)^2)/exp(N/2/Pi/s/m/log(N))+512/3*zeta(11/2)^2*Pi^7*m^(7/2)*sigma(N,0)/N^(11/2)/(1-o^2)^3,8/9*Pi^5*zeta(7/2)^2*m^(5/2)*s^2/N^(3/2)*log(N)^2);
}



\\ The bound for E^3 given by D-J in page 12, spliting via c<p^2*N, or not. The first sum is bounded by the integral,
\\ Take split=p^2.

E4(N,q,split,s)=
{if(s,,s=q^2/2/Pi);
16*Pi^3*(2/Pi*eulerphi(q)/N*(log(split)^2/2+log(N)*(log(split)+1))+s/6*log(N)/sqrt(N)*sigma(N,0)*(zeta(3/2)^2-sum(i=1,split,sigma(i,0)/sqrt(i^3))))
}


\\ Following D-J to get the lower bound for (a_1,L_\chi)_{p^2}.

F(p,q,s,split)=
{local(N);
N=p^2;
if(split,,split=N*p);
if(s,,s=q^2/2/Pi);
4*Pi*exp(-2*Pi/s/N/log(N))-E4(N,q,split,s)-E3(1,N,q,s)-E2(1,N,q,s)-E1(1,N,q,s)-bound1(1,N,q,s)
};

\\ Following D-J Lemma 7, to bound (a_1,L_\chi)_p.

F2(p,q,m)=2*sqrt(3)*sqrt(m)*sigma(m,0)/(1-exp(-2*Pi/q/sqrt(p)))*(4*Pi+16*zeta(3/2)^2*Pi^2/sqrt(p^3));

\\======================================================================
\\ This is the script to compute the bound. q = conductor K/Q. The first p where the output is positive is the right answer.

EllenbergBound(p,q,split)=F(p,q,q^2/2/Pi,split)-1/(p^2-1)*F2(p,q,p)-p/(p^2-1)*F2(p,q,1)






