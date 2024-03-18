function [C1,C2,t]=thin_lens_thickning_spherical_coma(K,ua,ua_,ub,ub_,ya,yb,n,S1_r,S2_r)

con=[1,3];
lb=[-500, 0.1];
up=[-500,20];
fun=@(con)thin_lens_thickning_spherical_solver(con,K,ua,ua_,ub,ub_,ya,yb,n,S1_r,S2_r);
ms = MultiStart('FunctionTolerance',1e-12,'UseParallel',true);
gs = GlobalSearch(ms);

problem = createOptimProblem('fmincon','x0',con,...
    'objective',fun,'lb',lb,'ub',ub);
x = run(gs,problem);
fprintf("Smallest Objective:%0.12e\n",fun(x));
C1=x(1);
t=x(2);
C2=(K-(n-1)*C1)/((1-n)+t/n*C1*(n-1)^2);



end


function deviation=thin_lens_thickning_spherical_solver(con,K,ua,ua_,ub,ub_,ya,yb,n,S1_r,S2_r)
%syms C1 t
C1=con(1);
t=con(2);
C2=(K-(n-1)*C1)/((1-n)+t/n*C1*(n-1)^2);
ya1=(ya+t*ua_)/(1-t*(1-n)*C2);
yb1=(yb+t*ub_)/(1-t*(1-n)*C2);
ua1=ua_+ya1*(1-n)*C2;
yb1=(yb+t*ub_)/(1-t*(1-n)*C2);
ub1=ub_+yb1*(1-n)*C2;%%%%%%%%%%%would this correct??????


A1=ua+ya*C1;
A_1=ub+yb*C1;
A2=ua_+ya1*C2;
A_2=ub_+yb1*C2;
S11=-A1^2*ya*(ua1/n-ua);
S12=-A2^2*ya1*(ua_/n-ua1);
S21=-A1*A_1*ya*(ua1/n-ua);
S22=-A2*A_2*ya1*(ua_/n-ua1);
S1=S11+S12;
S2=S21+S22;

w_expect=[S1_r,S2_r];
deviation=(([S1,S2]-(w_expect)));

%deviation(4)=deviation(4)*10;
deviation=deviation./(w_expect);
deviation=sum(abs(deviation));
end

