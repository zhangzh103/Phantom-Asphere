function result=singlet_aberration_contribution(C1,C2,ya,yb,ua,ua_,ub,n,phi,NA)%pezval sum is not correct
% m=ua/ua_;
% Y=(1+m)/(1-m);
% X=(C1+C2)/(C1-C2);
% H=1*ua*yb-1*ub*ya;
% A=(n+2)/(n*(n-1)^2);
% B=4*(n+1)/(n*(n-1));
% C=(3*n+2)/n;
% D=n^2/(n-1)^2;
% E=(n+1)/(n*(n-1));
% F=(2*n+1)/n;
% gamma1=A*X^2-B*X*Y+C*Y^2+D;
% gamma2=E*X-F*Y;
% w1=1/32*ya^4*phi^3*gamma1;
% w2=-1/4*H*ya^2*phi^2*gamma2;
% w3=1/2*H^2*phi;
% w4=1/4*H^2*phi/n;
% result=[w1,w2,w3,w4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0;
ua1=(ua-ya*(n-1)*C1)/n;
ya1=ya+n*ua1*t/n;
ua_=n*ua1-ya1*(1-n)*C2;
ub1=(ub-yb*(n-1)*C1)/n;
yb1=yb+n*ub1*t/n;
ub_=n*ub1-yb1*(1-n)*C2;
H=ua*yb-ub*ya;
A1=ua+ya*C1;
A_1=ub+yb*C1;
A2=ua_+ya1*C2;
A_2=ub_+yb1*C2;
P1=C1*(1/n-1);
P2=C2*(1-1/n);
P=[P1,P2];
inc_un1=(ua1/n-ua);
inc_un2=(ua_-ua1/n);
inc_un1_=(ub1/n-ub);
inc_un2_=(ub_-ub1/n);
inc_un_=[inc_un1_,inc_un2_];
A=[A1,A2];
A_=[A_1,A_2];
y=[ya,ya1];
y_=[yb,yb1];
inc_un=[inc_un1,inc_un2];

S1=-A.^2.*y.*inc_un;
S2=-A.*A_.*y.*inc_un;
S3=-A_.^2.*y.*inc_un;
S4=-H^2.*P;
S5=-A_./A.*(H^2.*P+A_.^2.*y.*inc_un);
S6=-A_.^2.*y_.*inc_un_;

S1=sum(S1);
S2=sum(S2);
S3=sum(S3);
S4=sum(S4);
S5=sum(S5);
S6=sum(S6);

W040=1/8*S1;
W131=1/2*S2;
W222=1/2*S3;
W220=1/4*(S4+S3);
w(4)=W220/NA*2;
w(1)=W040/NA*4;
w(2)=W131/NA*3;
w(3)=W222/NA*2+w(4);
w(4)=S4/4/NA*2;%why???


phi1=C1*(n-1);
phi2=C2*(1-n);
Phi=phi1+phi2-t*phi1*phi2;
delta=t*phi2/Phi;
delta_=-phi1*t*n/Phi;

result=[W040,W131,W222,S4/4,delta,delta_];
A=1;
end