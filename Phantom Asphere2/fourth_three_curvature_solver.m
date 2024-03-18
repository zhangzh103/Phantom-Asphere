function fourth_three_curvature_solver(C1s,C2s,Phi,dj,inc_n,ya,ya_,yb,ua,ua_,ub)
f1=-700;
f2=700;
f3=-700;
t1=0;
t2=0;
%C4=-700;
con=[f1 f2 f3 t1 t2];%initinal condition
%options=optimset('Display','iter', ...
%'PlotFcns',@optimplotfval);%Set coeffieicent
%options.FunctionTolerance = 1e-50;
fun=@(con)two_singlet_evaluater(C1s,C2s,con,Phi,dj,inc_n,ya,ya_,yb,ua,ua_,ub);
%bestRBF=fminsearch(fun,con,options);

%退火
%options.ReannealInterval = 100;
%options.ObjectiveLimit=0.001;
ub=[10000 10000 10000 10 10];
lb=[-10000 -10000 -10000 0 0];
%[bestRBF,fval,exitFlag,output]=simulannealbnd(fun,con,lb,ub,options)
%fprintf('The best function value found was : %g\n', fval);

%GS
ms = MultiStart('FunctionTolerance',1e-12,'UseParallel',true);
gs = GlobalSearch(ms)
problem = createOptimProblem('fmincon','x0',con,...
    'objective',fun,'lb',lb,'ub',ub);
x = run(gs,problem)


C1=1/x(1);
C2=1/x(2);
C3=1/x(3);
t1=x(4);
t2=x(5);
%C4=1/f_all(4);
phi1=0.5*C1-0.5*C2-t1*(-0.5*C2)*(0.5*C1);
%phi2=0.5*C3-0.5*C4;
phi2=Phi-phi1;
C4=(phi2-0.5*C3)/(-0.5+t2*0.5*C3);
f1=1/C1;
f2=1/C2;
f3=1/C3;
f4=1/C4;
FOCAL=[f1 f2 f3 f4];
fprintf("R1:%0.4f, R2:%0.4f, R3:%0.4f, R4:%0.4f, t1:%0.4f, t2:%0.4f\n",f1,f2,f3,f4,t1,t2)
fprintf("Smallest Objective:%0.12e\n",fun(x))
end