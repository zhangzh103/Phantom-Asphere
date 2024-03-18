function deviation=two_singlet_evaluater(C1s,C2s,C_all,Phi,dj,inc_n,ya,yb,ua,ua_,ub)
n=1.5;
C1=C_all(1);
C2=C_all(2);
C3=C_all(3);
%C4=C_all(4);
phi1=0.5*C1-0.5*C2;
phi2=Phi-0.5*C1+0.5*C2;
C4=C3-2*phi2;

ua2=ua-ya*phi1;
m1=ua/ua2;
m2=ua2/ua_;
Y1=(1+m1)/(1-m1);
Y2=(1+m2)/(1-m2);
X1=(C1+C2)/(C1-C2);
X2=(C3+C4)/(C3-C4);
H=1*ua*yb-1*ub*ya;
A=(n+2)/(n*(n-1)^2);
B=4*(n+1)/(n*(n-1));
C=(3*n+2)/n;
D=n^2/(n-1)^2;
E=(n+1)/(n*(n-1));
F=(2*n+1)/n;

if isnan(X1) || isinf(X1)
    X1=0;
end
if isnan(X2) || isinf(X2)
    X2=0;
end
if isnan(Y1) || isinf(Y1)
    Y1=0;
end
if isnan(Y2) || isinf(Y2)
    Y2=0;
end
gamma1_1=A*X1^2-B*X1*Y1+C*Y1^2+D;
gamma2_1=E*X1-F*Y1;

gamma1_2=A*X2^2-B*X2*Y2+C*Y2^2+D;
gamma2_2=E*X2-F*Y2;

w1_1=1/32*ya^4*phi1^3*gamma1_1;
w2_1=-1/4*H*ya^2*phi1^2*gamma2_1;
w3_1=1/2*H^2*phi1;
w4_1=1/4*H^2*phi1/n;

w1_2=1/32*ya^4*phi2^3*gamma1_2;
w2_2=-1/4*H*ya^2*phi2^2*gamma2_2;
w3_2=1/2*H^2*phi2;
w4_2=1/4*H^2*phi2/n;

w1=w1_1+w1_2;
w2=w2_1+w2_2;
w3=w3_1+w3_2;
w4=w4_1+w4_2;

w_expect=singlet_aberration_contribution(C1s,C2s,ya,yb,ua,ua_,ub,n,Phi)+asphere_aberration_contribution(dj,inc_n,ya,yb);
deviation=sum(abs([w1,w2,w3,w4]-w_expect).^0.5);
deviation=1/deviation*-1;
if isnan(deviation)
    fprintf("ERROR!!!")
    deviation=-1000000000000000000000;
end
end