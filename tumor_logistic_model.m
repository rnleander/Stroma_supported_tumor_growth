function dydt =tumor_logistic_model(t,y,lambda,K,kp)
%%%%%Adaptation of tumor model from Ribba et al
%%% P is the part of the tumor radius sourced from proliferative tissue
%%% Q is the part of the tumor radius sourced from quiescent tissue
%%% It is is assumed that tissue moves between the two compartments with
%%% rates kp and kq, and P increases logistically with per capita rate
%%% lambda and maximal radius K

P=y(1);
Q=y(2);
R=P+Q;

dPdt = lambda*P*(1-R/K)-kp*P;
dQdt = kp*P;

dydt=[dPdt;dQdt];

end