function dydt =tumor_growth_model_dC(t,y,mu,d,C,CG_ratio)
%%%%%The tumor growth model without immune system activation

DG=C/CG_ratio;
DC=C;

R=y;

if R>d
    dydt=d*mu*DG*(R^2-d*R+(1/3)*d^2)/(DC*R^2-2*d*(DC-DG)*R+(DC-DG)*d^2);
end

if R<=d
    dydt=mu*R/3;
end

end