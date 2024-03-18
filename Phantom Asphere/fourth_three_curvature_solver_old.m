function fourth_three_curvature_solver(Phi,dj,inc_n,ya,yb,ua,ua_,ub)
syms C1 C2
C4=0;
C3=0;
n=1.5;
%m=ua/ua_;
phi1=0.5*C1+0.5*C2;
phi2=0.5*C3+0.5*C4;
%s_F=(0.5*C1+0.5*C2)/(0.5*C3+0.5*C4);

%Y=(1+m)/(1-m);
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

gamma1_1=A*X1^2-B*X1*Y1+C*Y1^2+D;
gamma2_1=E*X1-F*Y1;

gamma1_2=A*X2^2-B*X2*Y2+C*Y2^2+D;
gamma2_2=E*X2-F*Y2;

w1_1=1/32*ya^4*phi1^3*gamma1_1;
w2_1=-1/4*H*ya^2*phi1^2*gamma2_1;
w3_1=1/2*H^2*phi1;
w4_1=1/4*H^2*phi1/n;

w1_2=1/32*ya^4*phi1^3*gamma1_2;
w2_2=-1/4*H*ya^2*phi1^2*gamma2_2;
w3_2=1/2*H^2*phi1;
w4_2=1/4*H^2*phi1/n;

w1=w1_1+w1_2;
w2=w2_1+w2_2;
w3=w3_1+w3_2;
w4=w4_1+w4_2;

w_expect=singlet_aberration_contribution(C1,C2,ya,yb,ua,ua_,ub,n,Phi)+asphere_aberration_contribution(dj,inc_n,ya,yb);
equation=[w1==w_expect(1),w2==w_expect(2),w3==w_expect(3),w4==w_expect(4)];
solve(equation)
end