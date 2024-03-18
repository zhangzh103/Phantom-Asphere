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
if isnan(sum(result)) || isinf(sum(result))
    X1=0;
    fprintf("Nan Present")
end

end