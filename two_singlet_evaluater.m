function deviation=two_singlet_evaluater(C1s,C2s,ft_all,Phi,dj,inc_n,ya,ya_,yb,ua,ua_,ub)
n=1.5;
C1=1/ft_all(1);
C2=1/ft_all(2);
C3=1/ft_all(3);
t1=ft_all(4);
t2=ft_all(5);
%C4=1/f_all(4);
phi1=0.5*C1-0.5*C2-t1*(-0.5*C2)*(0.5*C1);
%phi2=0.5*C3-0.5*C4;

%phi2=Phi-phi1;%cannot use this any more, use ua' and brick equations solve for C4 to make sure the Phi keep the same
%C4=(phi2-0.5*C3)/(-0.5+t2*0.5*C3);

%new
phi_s1=0.5*C1;
phi_s2=-0.5*C2;
phi_s3=0.5*C3;
ua2=ua-ya*phi_s1;
ya2=ya+ua2*t1;
ua3=ua2-ya2*phi_s2;
ya3=ya2;
ua4=ua3-ya3*phi_s3;
ya4=ya+ua3*t2;

C4=-2*(1.5*ua4-ua_)/ya4;
phi2=0.5*C3-0.5*C4-t2*(-0.5*C4)*(0.5*C3);


end