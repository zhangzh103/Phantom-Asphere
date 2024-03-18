function phantom_asphere_ver2()%%%%wronge!!! check the paraxial ray tracing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yb=0;
ya=1;
ub=0.0098333333333333333;
ub_=0.009506;
yb_=0.065556;
ya_=0.966667;
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
dj=0.001;
dj=+7.55053978932568e-07;
%dj=0
inc_n=0.5;
result=asphere_aberration_contribution(dj,inc_n,ya,yb);
[S1,S2,S3,S4,S5,S6,delta,delta_]=singlet_aberration_contribution(C1,C2,ya,yb,ua,ub,n,NA,t);
%[W040,W131,W222,W220]=S_2_W(S1,S2,S3,S4);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%


w=[0,0,0,0];
w_expect=[S1,S2,S3,S4]+result;
W040=1/8*w_expect(1);
W131=1/2*w_expect(2);
W222=1/2*w_expect(3);
W220=1/4*(w_expect(4)+w_expect(3));
w(4)=W220/NA*2;
w(1)=W040/NA*4;
w(2)=W131/NA*3;
w(3)=W222/NA*2+w(4);
w(4)=S4/4/NA*2;%why???
[C11,C21,t1,C12,C22,t2,d,deviationn]=thin_lens_thickning_spherical_coma_m2(ua,ua_,ub,ub_,ya,ya_,yb,yb_,n,w_expect,NA);
fprintf("R1:%0.4f, R2:%0.4f, R3:%0.4f, R4:%0.4f, t1:%0.4f, t2:%0.4f, d:%0.4f\n",1/C11,1/C21,1/C12,1/C22,t1,t2,d);
A=1;
end

function [C11,C21,t1,C12,C22,t2,d,deviation]=thin_lens_thickning_spherical_coma_m2(ua,ua_,ub,ub_,ya,ya_,yb,yb_,n,w_expect,NA)

con=[0.01,0.01,3, 0.01, 0.01, 3,1];%C11 C21 t1; C21 C22 t2
lb=[-5,-5,0.0001,-5,-5,0.0001,0.00001];
ubound=[5,5,10,5,5,10,10];

fun=@(con)thin_lens_thickning_spherical_solver_m2(con,ua,ua_,ub,ub_,ya,ya_,yb,yb_,n,w_expect,NA);
ms = MultiStart('FunctionTolerance',1e-12,'UseParallel',true);
gs = GlobalSearch(ms);

problem = createOptimProblem('fmincon','x0',con,...
    'objective',fun,'lb',lb,'ub',ubound);
x = run(gs,problem);
fprintf("Smallest Objective:%0.12e\n",fun(x));
deviation=fun(x);
C11=x(1);
C21=x(2);
t1=x(3);
C12=x(4);
C22=x(5);
t2=x(6);
d=x(7);





end


function deviation=thin_lens_thickning_spherical_solver_m2(con,ua,ua_,ub,ub_,ya,ya_,yb,yb_,n,w_expect,NA)%need one more!!! d 
C11=con(1);
t1=con(3);
C21=con(2);
d=con(7);
%C21=(K-(n-1)*C1)/((1-n)+t/n*C1*(n-1)^2);
ua1=(ua-ya*(n-1)*C11)/n;
ya1=ya+n*ua1*t1/n;
ua2=n*ua1-ya1*(1-n)*C21;
ub1=(ub-yb*(n-1)*C11)/n;
yb1=yb+n*ub1*t1/n;
ub2=n*ub1-yb1*(1-n)*C21;

ya2=ya1+ua2*d;
yb2=yb1+ub2*d;
C12=con(4);
t2=con(6);
C22=con(5);
ua3=(ua2-ya2*(n-1)*C12)/n;
ub3=(ub2-yb2*(n-1)*C12)/n;
%ub3=(n*ub2-yb2*(1-n)*C12)/n;
ya3=ya2+ua3*t2;
yb3=yb2+ub3*t2;

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
height_diff=[(ya3-ya_),(yb3-yb_)];
height_diff=height_diff/[ya_,yb_];
height_diff=sum(abs(height_diff));
%deviation(4)=deviation(4)*10;
phi1=(n-1)*C11;
phi2=(1-n)*C21;
%Phi1=phi1+phi2-t*phi1*phi2/n;

deviation=deviation./(w_expect);
deviation=sum(abs(deviation))+abs(angle_dif)+height_diff;
%deviation=sum(abs(deviation));
end

function [S1,S2,S3,S4,S5,S6,delta,delta_]=singlet_aberration_contribution(C1,C2,ya,yb,ua,ub,n,NA,t)%pezval sum is not correct
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

function result=asphere_aberration_contribution(dj,inc_n,ya,yb)
a=8*dj*ya^4*inc_n;
ybar_y=yb/ya;
S1=a;
S2=ybar_y*a;
S3=ybar_y^2*a;
S4=0;
result=[S1,S2,S3,S4];

end