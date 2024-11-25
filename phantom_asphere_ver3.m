
function phantom_asphere_ver3()%%%%%%%%%%%%%%one more thing!!! principle plane modification
addpath("Codev_Dependent")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

surface_num=11;%surface number of the asphere
[R1,R2,t,n,ya,ya_,yb,yb_,ua,ua_,ub,ub_,dj,NA]=basic_data_generator(surface_num);%grab data
C1=1/R1;
C2=1/R2;
dj=dj/100000000;
inc_n=n-1;

result=asphere_aberration_contribution(dj,inc_n,ya,yb);%calculate the original singlet contribution
[S1,S2,S3,S4,S5,S6,delta,delta_]=singlet_aberration_contribution(C1,C2,ya,yb,ua,ub,n,NA,t);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%


w=[0,0,0,0];
w_expect=[S1,S2,S3,S4]+result;%pure singlet contribution
W040=1/8*w_expect(1);
W131=1/2*w_expect(2);
W222=1/2*w_expect(3);
W220=1/4*(w_expect(4)+w_expect(3));
w(4)=W220/NA*2;
w(1)=W040/NA*4;
w(2)=W131/NA*3;
w(3)=W222/NA*2+w(4);
w(4)=S4/4/NA*2;%why??? Not fully the same as the equation but correct
deviation=1;
while deviation>0.01%allowed error
[C11,C21,t1,C12,C22,t2,d,deviation,delta2,delta2_]=thin_lens_thickning_spherical_coma_m2(ua,ua_,ub,ub_,ya,ya_,yb,yb_,n,w_expect,NA);
end
%The lens go like this way:
%C11 t1 C21 d C12 t2 C22, the last air space comes fron the previosu
fprintf("R1:%0.4f, R2:%0.4f, R3:%0.4f, R4:%0.4f, t1:%0.4f, t2:%0.4f, d:%0.4f\n",1/C11,1/C21,1/C12,1/C22,t1,t2,d);
A=1;
end

function [C11,C21,t1,C12,C22,t2,d,deviation,delta,delta_]=thin_lens_thickning_spherical_coma_m2(ua,ua_,ub,ub_,ya,ya_,yb,yb_,n,w_expect,NA)

con=[0.01,0.01,3, 0.01, 0.01, 3,1];%C11 C21 t1; C21 C22 t2 d
lb=[-20,-20,0.0001,-20,-20,0.0001,0.00001];
ubound=[20,20,10,20,20,10,10];
%%%%%%%%%%%%%%%%%%%%%%%%R
%con=[100,-100,3, 100, -100, 3,1];%C11 C21 t1; C21 C22 t2
%lb=[-2000,-2000,0.0001,-2000,-2000,0.0001,0.00001];
%ubound=[2000,2000,10,2000,2000,10,10];
%%%%%%%%%%%%%

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
[ua1,ua2,ua3,ua4,ub1,ub2,ub3,ub4,ya1,ya2,ya3,yb1,yb2,yb3]=paraxial_ray_trace_two_thick(C11,C21,C12,C22,t1,t2,d,0,ub,ya,yb,n);
delta=(ya3-ya)/tan(ua4);
[ua1,ua2,ua3,ua4,ub1,ub2,ub3,ub4,ya1,ya2,ya3,yb1,yb2,yb3]=paraxial_ray_trace_two_thick(-C22,-C12,-C21,-C11,t2,t1,d,0,ub,ya,yb,n);
delta_=-(ya3-ya)/tan(ua4);




end


function deviation=thin_lens_thickning_spherical_solver_m2(con,ua,ua_,ub,ub_,ya,ya_,yb,yb_,n,w_expect,NA)%need one more!!! d 
C11=con(1);
t1=con(3);
C21=con(2);
d=con(7);
C12=con(4);
t2=con(6);
C22=con(5);
[ua1,ua2,ua3,ua4,ub1,ub2,ub3,ub4,ya1,ya2,ya3,yb1,yb2,yb3]=paraxial_ray_trace_two_thick(C11,C21,C12,C22,t1,t2,d,ua,ub,ya,yb,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%from singlet_aberration_contribution
%Use the ray trace data to calculate aberration contribution
[S1_1,S2_1,S3_1,S4_1,S5_1,S6_1,delta,delta_]=singlet_aberration_contribution(C11,C21,ya,yb,ua,ub,n,NA,t1);%%%%%%%%%%%%%%%%%%%%%%%
[S1_2,S2_2,S3_2,S4_2,S5_2,S6_2,delta,delta_]=singlet_aberration_contribution(C12,C22,ya2,yb2,ua2,ub2,n,NA,t2);
w_real=[S1_1+S1_2,S2_1+S2_2,S3_1+S3_2,S4_1+S4_2];%aberration contribtion of the two lenses
deviation=((w_real-(w_expect)));
angle_dif=[ua4-ua_,ub4-ub_];%The angle difference compare with the asphere
angle_dif=angle_dif/[ua_,ub_];
height_diff=[(ya3-ya_),(yb3-yb_)];
height_diff=height_diff/[ya_,yb_];
height_diff=sum(abs(height_diff));
phi1=(n-1)*C11;
phi2=(1-n)*C21;

deviation=deviation./(w_expect);
deviation=sum(abs(deviation))+abs(angle_dif)+height_diff;
end

function [S1,S2,S3,S4,S5,S6,delta,delta_]=singlet_aberration_contribution(C1,C2,ya,yb,ua,ub,n,NA,t)%pezval sum is not correct
%calculate the aberration contribution from one siglet
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
delta=t*phi2/n/Phi;
delta_=-phi1*t/n/Phi;

result=[S1,S2,S3,S4,S5,S6,delta,delta_];
end
function [ua1,ua2,ua3,ua4,ub1,ub2,ub3,ub4,ya1,ya2,ya3,yb1,yb2,yb3]=paraxial_ray_trace_two_thick(C11,C21,C12,C22,t1,t2,d,ua,ub,ya,yb,n)
%paraxial ray tracing for two thick lens with airspace
ua1=(ua-ya*(n-1)*C11)/n;
ya1=ya+n*ua1*t1/n;
ua2=n*ua1-ya1*(1-n)*C21;
ub1=(ub-yb*(n-1)*C11)/n;
yb1=yb+n*ub1*t1/n;
ub2=n*ub1-yb1*(1-n)*C21;
ya2=ya1+ua2*d;%air space
yb2=yb1+ub2*d;
ua3=(ua2-ya2*(n-1)*C12)/n;
ub3=(ub2-yb2*(n-1)*C12)/n;
ya3=ya2+ua3*t2;
yb3=yb2+ub3*t2;
ua4=n*ua3-ya3*(1-n)*C22;
ub4=n*ub3-yb3*(1-n)*C22;

end
function result=asphere_aberration_contribution(dj,inc_n,ya,yb)
%do not use any more, grab the data directly from the Codev
a=8*dj*ya^4*inc_n;
ybar_y=yb/ya;
S1=a;
S2=ybar_y*a;
S3=ybar_y^2*a;
S4=0;
result=[S1,S2,S3,S4];

end