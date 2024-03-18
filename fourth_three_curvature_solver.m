function fourth_three_curvature_solver(C1s,C2s,Phi,dj,inc_n,ya,yb,ua,ua_,ub)
C1=-700;
C2=700;
C3=-700;
%C4=-700;
con=[C1 C2 C3];%initinal condition
%options=optimset('Display','iter', ...
%'PlotFcns',@optimplotfval);%Set coeffieicent
%options.FunctionTolerance = 1e-50;
fun=@(con)two_singlet_evaluater(C1s,C2s,con,Phi,dj,inc_n,ya,yb,ua,ua_,ub);
%bestRBF=fminsearch(fun,con,options);

%退火
%options.ReannealInterval = 100;
%options.ObjectiveLimit=0.001;
ub=[10000 10000 10000];
lb=[-10000 -10000 -10000];
%[bestRBF,fval,exitFlag,output]=simulannealbnd(fun,con,lb,ub,options)
%fprintf('The best function value found was : %g\n', fval);

%GS
ms = MultiStart('FunctionTolerance',1e-12,'UseParallel',true);
gs = GlobalSearch(ms)
problem = createOptimProblem('fmincon','x0',con,...
    'objective',fun,'lb',lb,'ub',ub);
x = run(gs,problem)


f1=x(1);
f2=x(2);
f3=x(3);
%f4=x(4);
C3=1/f3;
phi2=Phi-0.5*1/f1+0.5*1/f2;
C4=C3-2*phi2;
f4=1/C4;
FOCAL=[f1 f2 f3 f4];
fprintf("R1:%0.4f, R2:%0.4f, R3:%0.4f, R4:%0.4f\n",f1,f2,f3,f4)
fprintf("Smallest Objective:%0.12e\n",fun(x))
end