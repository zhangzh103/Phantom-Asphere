function main()

delta=1.64972;
delta_=-7.98994;
d=12.35;
d1=d/2;
d2=d-d1;
K=-0.107904;
ua=-0.0344412;
ya=1.46209;
ub=0.893084;
yb=-8.07159;
S1=-0.000957597;
S2=-0.00249863;
S3=-0.01697;
S4=-0.0649725;
S5=0.825262;
S6=-3.33372;

% ua=-0.161533;
% ya=2.98998;
% ub=0.338591;
% yb=-2.1769;
% S1=-0.00966773;
% S2=0.0117646;
% S3=-0.0381275;
% S4=-0.0268528;
% S5=0.0777384;
% S6=-0.0504529;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yb=0;
ya=1;
ub=0.0098333333333333333;
ua=0;
ua_=-0.010000;
n=1.5;
C1=1/100;
C2=1/-100;
%phi=1/100;
NA=0.00983333333333;
t=10;
phi1=C1*(n-1);
phi2=C2*(1-n);
d=10;
d1=2.5;
K=phi1+phi2-t*phi1*phi2/n;
[S1,S2,S3,S4,S5,S6,delta,delta_]=singlet_aberration_contribution(C1,C2,ya,yb,ua,ub,n,NA,t);
[W040,W131,W222,W220]=S_2_W(S1,S2,S3,S4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%



K_all=thin_lens_power_calculatora(delta,delta_,K,d1,d);
[S,n,ya,yb,ua,ub]=third_order_aberration_calculator(ua,ub,ya,yb,K_all,d1,d2,S1,S2,S3,S4,S5,S6);
S1_1=S(1);
S2C_1=S(2);
S1_2=S(3);
S2C_2=S(4);
S1_3=S(5);
S2C_3=0;
ya1=ya(2);
ya2=ya(3);
ya=ya(1);
yb1=yb(2);
yb2=yb(3);
yb=yb(1);

ub1=ub(2);
ub2=ub(3);
ub3=ub(4);
ub=ub(1);

ua1=ua(2);
ua2=ua(3);
ua3=ua(4);
ua=ua(1);
%[C13,C23]=turn_real_lens(K_all(3),S1_3,S2C_3+1e-12,n,ua2,ua3,ub2,ub3,ya2,yb2);%don't know why cannot have coma as zero
%[C11,C21]=turn_real_lens(K_all(1),S1_1,S2C_1,n,ua,ua1,ub,ub1,ya,yb);
%[C12,C22]=turn_real_lens(K_all(2),S1_2,S2C_2,n,ua1,ua2,ub1,ub2,ya1,yb1);
d2=d-d1;
K_cal=K_all(1)+K_all(2)+K_all(3)-d1*K_all(1)*K_all(2)-d1*K_all(1)*K_all(3)-d2*K_all(1)*K_all(3)-d2*K_all(2)*K_all(3)+d1*d2*K_all(1)*K_all(2)*K_all(3);
deviation=1;
while deviation>0.2
[C11,C21,t1,deviation]=thin_lens_thickning_spherical_coma(K_all(1),ua,ua1,ub,ub1,ya,yb,n,S1_1,S2C_1,NA);
end
%fprintf("n:%0.4f; R1:%0.4f; R2:%0.4f; R3:%0.4f; R4:%0.4f; R5:%0.4f; R6:%0.4f\n",n,1/C11,1/C21,1/C12 ...
%    ,1/C22,1/C13,1/C23)
phi1=C11*(n-1);
phi2=C21*(1-n);
Phi1=phi1+phi2-t1*phi1*phi2/n;
deviation=1;
ua1=(ua-ya*(n-1)*C11)/n;
ya1=ya+n*ua1*t/n;
ua_=n*ua1-ya1*(1-n)*C21;
ub1=(ub-yb*(n-1)*C11)/n;
yb1=yb+n*ub1*t/n;
ub_=n*ub1-yb1*(1-n)*C21;
H=ua*yb-ub*ya;
singlet_aberration_contribution(C11,C21,ya,yb,ua,ub,n,NA,t);
while deviation>0.2
[C12,C22,t2,deviation]=thin_lens_thickning_spherical_coma(K_all(2)+K_all(3)+K_all(3)*K_all(2),ua_,ua3,ub_,ub3,ya1,yb1,n,S1_2+S1_3,S2C_2+S2C_3,NA);
end
phi1=C12*(n-1);
phi2=C22*(1-n);
Phi2=phi1+phi2-t2*phi1*phi2/n;
A=1;
end

function result=thin_lens_power_calculatora(delta,delta_,K,d1,d)
%delta front principle plane p
%delta_ back principle plane p'
%d1: distance between first two element
%d2: distance between second two element
%K is the power of the system
%K1 is the power of the first thin lens
d2=d-d1;
K1=((delta-d1)*delta_*K-(d-delta+delta_))/((delta-delta_+delta*delta_*K)*d1);
K2=(d-(delta-delta_+delta*delta_*K))/(d1*d2);
K3=((d2+delta_)*delta*K-(d-delta+delta_))/((delta-delta_+delta*delta_*K)*d2);
result=[K1,K2,K3];
end

function [S,n,ya,yb,ua,ub]=third_order_aberration_calculator(ua,ub,ya,yb,K_all,d1,d2,S1,S2,S3,S4,S5,S6)
%h_bar: chief ray
%h: marginal ray
n=1.5;
K1=K_all(1);
K2=K_all(2);
K3=K_all(3);
%%first lens
ua1=ua-ya*K1;
ub1=ub-yb*K1;
ya1=ya+ua1*d1;
yb1=yb+ub1*d1;
ua2=ua1-ya1*K2;
ub2=ub1-yb1*K2;
ya2=ya1+ua2*d2;
yb2=yb1+ub2*d2;
ua3=ua2-ya2*K3;
ub3=ub2-yb2*K3;
E1=yb/ya;
E2=yb1/ya1;
E3=yb2/ya2;
H=ua*yb-ub*ya;
E=[E1 E2 E3];
S7=[H*(ua1^2-ua^2) H*(ua2^2-ua1^2) H*(ua3^2-ua2^2)];
S8=[H*(ua1*ub1-ua*ub) H*(ua2*ub2-ua1*ub1) H*(ua3*ub3-ua2*ub2)];
S9=[H*(ub1^2-ub^2) H*(ub2^2-ub1^2) H*(ub3^2-ub2^2)];
%S4_sudo=[H^2*K1 H^2*K2 H^2*K3];
%n=sum(S4_sudo)./S4;
n=1.5;
S4=[K1/(K1+K2+K3)*S4,K2/(K1+K2+K3)*S4,K3/(K1+K2+K3)*S4 ];
b1=S1;
b2=S2;
b3=S3-H^2*sum(K_all);
b4=S5-sum(E.*(3*H^2.*K_all+S4));
b5=S6-sum(6*H^2.*E.^2.*K_all+2*E.^2.*S4-E.^3.*S7+3*E.^2.*S8-3*E.*S9);

a16=(2*(2*E1-E2-E3)*(E2-E3)^2)/((E1-E2)^3*(E1-E3));
a26=-(E2-E3)^2/(E1-E2)^2;
a36=(2*(E1^2-E3^2)*(E1-2*E2-E3)+8*E1*E2*E3)/((E1-E2)^3*(E2-E3));
a36=(2*(E1^2-E3^2)*(E1-2*E2-E3)+8*E1*E2*E3)/((E1-E2)^3*(E2-E3));
a36=-1.55293;
a46=-(E1-E3)^2/(E1-E2)^2;
a56=-2*(E1+E2-2*E3)/((E1-E3)*(E2-E3));
g1=1/((E1-E2)^3*(E1-E3)^2)*(E2^2*E3*(2*E1*(-2*E1+E2) ...
    +E3*(3*E1-E2))*b1+2*E1*E2*((2*E1-E2)*(E2+2*E3) ...
    -3*E3^2)*b2+(-4*E1^2*(2*E2+E3)+(E1+E2) ...
    *(E2^2+2*E2*E3+3*E3^2))*b3 ...
    +2*(2*E1*(E1+E2)-(E2+E3)^2)*b4 ...
    +(-3*E1+E2+2*E3)*b5);

g2=1/((E1-E2)^2*(E1-E3))*(E1*E2^2*E3*b1-E2*(E2*E3 ...
    +E1*(E2+2*E3))*b2+(E1*(2*E2+E3)+E2*(E2+2*E3))*b3 ...
    -(E1+2*E2+E3)*b4+b5);
g3=1/((E1-E2)^3*(E2-E3)^2)*(E1^2*E3*(E2*(4*E2-3*E3) ...
    +E1*(-2*E2+E3))*b1+2*E1*E2*(E1^2-2*E1*E2+2*E1*E3 ...
    -4*E2*E3+3*E3^2)*b2+(-E1^3+E1*E2*(-E1+8*E2 ...
    -2*E3)+E1*E3*(-2*E1-3*E3)+E2*E3*(4*E2-3*E3))*b3 ...
    +2*(E1^2-2*E1*E2-2*E2^2+2*E1*E3+E3^2)*b4 ...
    +(-E1+3*E2-2*E3)*b5);
g4=1/((E1-E2)^2*(E2-E3))*(E1^2*E2*E3*b1-E1*(2*E2*E3 ...
    +E1*(E2+E3))*b2+(E1^2+E2*E3+2*E1*(E2+E3))*b3 ...
    -(2*E1+E2+E3)*b4+b5);
g5=1/((E1-E3)^2*(E2-E3)^2)*(E1^2*E2^2*b1-2*E1*E2*(E1+E2)*b2 ...
    +(E1^2+4*E1*E2+E2^2)*b3-2*(E1+E2)*b4+b5);
A=[1 0 0 0 0 a16; ...
   0 1 0 0 0 a26; ...
   0 0 1 0 0 a36; ...
   0 0 0 1 0 a46; ...
   0 0 0 0 1 a56];
G=[g1;g2;g3;g4;g5];
S=pinv(A)*G;
%%S= S1_1, S2c_1, S1_2, S2C_2, S1_3, S2C_3
S=G;%%assume S2,3=0!!!
%S=G*inv(A);
ua=[ua ua1 ua2 ua3];
ub=[ub ub1 ub2 ub3];
ya=[ya ya1 ya2];
yb=[yb yb1 yb2];
%result=[S,n,ya,yb];

end

function [C1,C2,t,deviation]=thin_lens_thickning_spherical_coma(K,ua,ua_,ub,ub_,ya,yb,n,S1_r,S2_r,NA)

con=[0.01,0.01,3];
lb=[-5,-5,0.0001];
ubound=[5,5,20];

fun=@(con)thin_lens_thickning_spherical_solver(con,K,ua,ua_,ub,ub_,ya,yb,n,S1_r,S2_r,NA);
ms = MultiStart('FunctionTolerance',1e-12,'UseParallel',true);
gs = GlobalSearch(ms);

problem = createOptimProblem('fmincon','x0',con,...
    'objective',fun,'lb',lb,'ub',ubound);
%%%%%%%%%%%%%%%%%%%%%%analyser
% figure(1)
% clf
% C1=linspace(1,-1,1000);
% t=linspace(0.1,20,1000);
% result=zeros(1000,1000);
% for i=1:1000
%     for j=1:1000
%     result(i,j)=fun([C1(i),t(j)]);
%     end
% end
% [C1,t]=meshgrid(C1,t);
% mesh(C1,t,result)
% xlabel("C1");
% ylabel("t");
% zlabel("Deviation")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = run(gs,problem);
fprintf("Smallest Objective:%0.12e\n",fun(x));
deviation=fun(x);
C1=x(1);
t=x(3);
C2=(K-(n-1)*C1)/((1-n)+t/n*C1*(n-1)^2);



end


function deviation=thin_lens_thickning_spherical_solver(con,K,ua,ua_,ub,ub_,ya,yb,n,S1_r,S2_r,NA)
%syms C1 t
C1=con(1);
t=con(3);
%C2=con(2);
C2=(K-(n-1)*C1)/((1-n)+t/n*C1*(n-1)^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%from singlet_aberration_contribution

[S1,S2,S3,S4,S5,S6,delta,delta_]=singlet_aberration_contribution(C1,C2,ya,yb,ua,ub,n,NA,t);%%%%%%%%%%%%%%%%%%%%%%%

w_expect=[S1_r,S2_r];
deviation=(([S1,S2]-(w_expect)));

%deviation(4)=deviation(4)*10;
phi1=(n-1)*C1;
phi2=(1-n)*C2;
Phi=phi1+phi2-t*phi1*phi2/n;
deviation=deviation./(w_expect);
deviation=sum(abs(deviation))+abs((Phi-K)/K);
%deviation=sum(abs(deviation));
end



function [S1,S2,S3,S4,S5,S6,delta,delta_]=singlet_aberration_contribution(C1,C2,ya,yb,ua,ub,n,NA,t)%pezval sum is not correct
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
%t=0;
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
Phi=phi1+phi2-t*phi1*phi2/n;
delta=t*phi2/Phi;
delta_=-phi1*t*n/Phi;

result=[S1,S2,S3,S4,S5,S6,delta,delta_];
end
function [W040,W131,W222,W220]=S_2_W(S1,S2,S3,S4)
W040=1/8*S1;
W131=1/2*S2;
W222=1/2*S3;
W220=1/4*(S4+S3);

end

function result=asphere_aberration_contribution(dj,inc_n,ya,yb)
% w_1=dj.*inc_n.*ya.^4;
% w_2=-0.5*ya/yb*(8*dj*ya.^4*inc_n);
% w_3=-0.5*(ya./yb).^2*(8*dj*ya.^4*inc_n);
% w_4=-0.25*(ya./yb).^2*(8*dj*ya.^4*inc_n);

kj_bar=-4*dj*inc_n;
w_1=0.25*kj_bar*ya^4;
w_2=1/3*kj_bar*ya^3*yb;
w_3=1/3*kj_bar*ya^2*yb^2;
w_4=0;
result=[w_1,w_2,w_3,w_4];

end

function [C1,C2,t,deviation]=thin_lens_thickning_spherical_coma_m2(ua,ua_,ub,ub_,ya,yb,n,w_expect,NA)

con=[0.01,0.01,3, 0.01, 0.01, 3];%C11 C21 t1; C21 C22 t2
lb=[-5,-5,0.0001,-5,-5,0.0001];
ubound=[5,5,20,5,5,20];

fun=@(con)thin_lens_thickning_spherical_solver_m2(con,ua,ua_,ub,ub_,ya,yb,n,w_expect,NA);
ms = MultiStart('FunctionTolerance',1e-12,'UseParallel',true);
gs = GlobalSearch(ms);

problem = createOptimProblem('fmincon','x0',con,...
    'objective',fun,'lb',lb,'ub',ubound);
%%%%%%%%%%%%%%%%%%%%%%analyser
% figure(1)
% clf
% C1=linspace(1,-1,1000);
% t=linspace(0.1,20,1000);
% result=zeros(1000,1000);
% for i=1:1000
%     for j=1:1000
%     result(i,j)=fun([C1(i),t(j)]);
%     end
% end
% [C1,t]=meshgrid(C1,t);
% mesh(C1,t,result)
% xlabel("C1");
% ylabel("t");
% zlabel("Deviation")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = run(gs,problem);
fprintf("Smallest Objective:%0.12e\n",fun(x));
deviation=fun(x);
C1=x(1);
t=x(3);
C2=(K-(n-1)*C1)/((1-n)+t/n*C1*(n-1)^2);



end


function deviation=thin_lens_thickning_spherical_solver_m2(con,ua,ua_,ub,ub_,ya,yb,n,w_expect,NA)%need one more!!! d 
%syms C1 t
C11=con(1);
t1=con(3);
C21=con(2);
%C21=(K-(n-1)*C1)/((1-n)+t/n*C1*(n-1)^2);
ua1=(ua-ya*(n-1)*C11)/n;
ya1=ya+n*ua1*t/n;
ua2=n*ua1-ya1*(1-n)*C21;
ub1=(ub-yb*(n-1)*C11)/n;
yb1=yb+n*ub1*t/n;
ub2=n*ub1-yb1*(1-n)*C21;

ya2=ya+ua2*d;
yb2=yb+ub2*d;
C12=con(4);
t2=con(6);
C22=con(5);
d=con(7);
ua3=(n*ua2-ya2*(1-n)*C12)/n;
ub3=(n*ub2-yb2*(1-n)*C12)/n;
ya3=ya2+ua3*d;
yb3=yb2+ub3*d;

ua4=n*ua3-ya3*(1-n)*C22;
ub4=n*ub3-yb3*(1-n)*C22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%from singlet_aberration_contribution

[S1_1,S2_1,S3_1,S4_1,S5_1,S6_1,delta,delta_]=singlet_aberration_contribution(C11,C21,ya,yb,ua,ub,n,NA,t1);%%%%%%%%%%%%%%%%%%%%%%%
[S1_2,S2_2,S3_2,S4_2,S5_2,S6_2,delta,delta_]=singlet_aberration_contribution(C12,C22,ya2,yb2,ua2,ub2,n,NA,t2);
%w_expect=[S1_r,S2_r];
w_real=[S1_1+S1_2,S2_1+S2_2,S3_1+S3_2,S4_1+S4_2];
deviation=((w_real-(w_expect)));
angle_dif=[ua4-ua_,ub4-ub_];
angle_dif=angle_dif/[ua_,ub_];
%deviation(4)=deviation(4)*10;
phi1=(n-1)*C11;
phi2=(1-n)*C21;
Phi1=phi1+phi2-t*phi1*phi2/n;

deviation=deviation./(w_expect);
deviation=sum(abs(deviation))+abs(angle_dif);
%deviation=sum(abs(deviation));
end
