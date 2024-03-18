
%fun=@(con)two_singlet_evaluater(C1s,C2s,con,Phi,dj,inc_n,ya,yb,ua,ua_,ub);
%bestRBF=fminsearch(fun,con,options);
n=1.5;
syms C1 C2
C3=0;
phi1=0.5*C1-0.5*C2;
phi2=Phi-0.5*C1+0.5*C2;
C4=C3-2*phi2;

ua2=ua-ya*phi1;
m1=ua/ua2;
m2=ua2/ua_;

yb=0;
ya=1;
ub=0.04;
ua=0.002;
ua_=-0.008;
n=1.5;
C1=1/100;
C2=1/-100;
phi=1/100;
dj=0;
inc_n=0.5;
Phi=1/100;

range=linspace(20,400,1000);
f1=range;
%f2=range;
%[f1,f2]=meshgrid(range);
f2=-range;
result=zeros(1000,1000);
for i=1:1000
for j=1:1000
    result(i,j)=two_singlet_evaluater(C1,C2,[f1(i),f2(j)],Phi,dj,inc_n,ya,yb,ua,ua_,ub);
end
end
figure(1)
clf
[f1,f2]=meshgrid(f1,f2);
hold on
%plot3(f1,f2,normalize(result))
mesh(f1,f2,result)
xlabel("f1")
ylabel("f2")
hold off