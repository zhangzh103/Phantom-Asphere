%%%%%%%%Reference: Modelling the first- and third-order properties of conceptual thick lenses by using three air-spaced thin lenses
%%%%%can you generate a singlet with arbitary sperhical aberration and coma
function main_three_thin()
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

[C11,C21]=thin_lens_thickning_spherical_coma(K_all(1),ua,ua1,ub,ub1,ya,yb,n,S1_1,S2C_1);
%fprintf("n:%0.4f; R1:%0.4f; R2:%0.4f; R3:%0.4f; R4:%0.4f; R5:%0.4f; R6:%0.4f\n",n,1/C11,1/C21,1/C12 ...
%    ,1/C22,1/C13,1/C23)

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
n=1;
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
S4_sudo=[H^2*K1 H^2*K2 H^2*K3];
n=sum(S4_sudo)./S4;
S4=S4_sudo./n;
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
function [C1, C2]=turn_real_lens(K,S1,S2C,n,ua,ua_,ub,ub_,ya,yb)
% syms C1 C2
% phi=(n-1)*C1+(1-n)*C2;
% H=ua*yb-ub*ya;
% Y=(ua_+ua)/(ua_-ua);
% X=(C1+C2)/(C1-C2);
% A=(n+2)/(n*(n-1)^2);
% B=4*(n+1)/(n*(n-1));
% C=(3*n+2)/n;
% D=n^2/(n-1)^2;
% E=(n+1)/(n*(n-1));
% F=(2*n+1)/n;
% gamma1=A*X^2-B*X*Y+C*Y^2+D;
% gamma2=E*X-F*Y;
% w1=1/32*ya^4*K^3*gamma1;
% w2=-1/4*H*ya^2*K^2*gamma2;
% equ=[w1==S1, phi==K,w2==S2C];
% result=solve(equ);
C1=K/(n-1);
fun=@(C1)singlet_evaluater(C1,K,n,S1,S2C,ya,yb,ua,ua_,ub);
lb=[-500];
ub=[500];
ms = MultiStart('FunctionTolerance',1e-12,'UseParallel',true);
gs = GlobalSearch(ms);

problem = createOptimProblem('fmincon','x0',C1,...
    'objective',fun,'lb',lb,'ub',ub);
x = run(gs,problem);
C1=x;
C2=(K-(n-1))*x/(1-n);
fprintf("Smallest Objective:%0.12e\n",fun(x))
end
function deviation=singlet_evaluater(C1,K,n,S1,S2C,ya,yb,ua,ua_,ub)
C2=(K-(n-1)*C1)/(1-n);


m1=ua/ua_;

Y1=(1+m1)/(1-m1);

X1=(C1+C2)/(C1-C2);

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
gamma1_1=A*X1^2-B*X1*Y1+C*Y1^2+D;
gamma2_1=E*X1-F*Y1;



w1_1=1/32*ya^4*K^3*gamma1_1;
w2_1=-1/4*H*ya^2*K^2*gamma2_1;
w3_1=1/2*H^2*K;
w4_1=1/4*H^2*K/n;



w1=w1_1;
w2=w2_1;
w3=w3_1;
w4=w4_1;

w_expect=[S1,S2C];
deviation=(abs([w1,w2]-w_expect));

%deviation(4)=deviation(4)*10;
deviation=deviation./abs(w_expect);
deviation=sum(deviation);
%deviation=1/deviation*-1;
if isnan(deviation)
    fprintf("ERROR!!!")
    deviation=100;
end
end
