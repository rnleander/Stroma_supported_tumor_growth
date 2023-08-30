function dydt =tumor_allometric_growth_model_d_Rbeta(t,y,mu,rho,CG_ratio,beta)
%%%%%The tumor growth model without immune system activation

R=y;
num=3*rho*(R.^(2+beta))+3*(rho^2).*(R^(2*beta+1))+(rho^3)*(R^(3*beta));
denom=3*rho*(2+beta)*(R.^(1+beta))+3*(2*beta+1)*(rho^2)*(R.^(2*beta))+3*beta*(rho^3)*R.^(3*beta-1)+3*CG_ratio*R^2;
g=num/denom;
dydt=mu*g;
end