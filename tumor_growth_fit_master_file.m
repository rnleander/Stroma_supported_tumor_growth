%Fit tumor growth data to all models, Master file%
function[mu,d,C,CG_ratio_d_constant,R0_d_C, mu_d_Rbeta,rho_d_Rbeta,CG_ratio_d_Rbeta,R0_beta,a,b,V0_exp]=tumor_growth_fit_master_file(mu_d_constant0,d00,C00,CG_ratio_d_constant00,R0_d_C00,mu_d_Rbeta0,rho_d_Rbeta0,CG_ratio_d_Rbeta0,R0_beta00,b00,V0_exp00)


%select error type
%error_type='rel_error';
error_type='meas_error';
%error_type='abs_error';
%error_type='LL_error';

%load tumor data to produce a matrices of observation times and volumes%

%read_mouse_data
%data_source='JIMT-1';

%read_glioma_data
%data_source='glioma';

%read_lung_cancer_data
%data_source='NSC lung cancer';

read_hamartoma_data
data_source='hamartoma';




nP=size(TV);
%nP number of patients%
nP=nP(2);


%find the data that are suitable for analysis
j=1;
for i=1:nP
    %remove NaN entries or entries were V==0
   
    nan_index = isnan(TV(:,i));
    valid_index=1-nan_index;
    nonzero_index=TV(:,i)~=0;
    valid_index=valid_index.*nonzero_index;

    valid_index=logical(valid_index);

    pdv=TV(valid_index,i);
    pdt=TT(valid_index,i);

    n=length(pdt);
%only keep data series with at least 5 observations
    if n>=5
    data_time.(sprintf('t%d',j))=pdt(1:n);
    data_v.(sprintf('v%d',j))=pdv(1:n);
    j=j+1;
    end
end
%find the number of suitable data sets;
nP=j-1;

%save error of fit to each patient
error_d_constant= zeros(1,nP);
error_exponential= zeros(1,nP);
error_d_R2=zeros(1,nP);
error_d_R=zeros(1,nP);
error_d_R3=zeros(1,nP);
error_d_Rbeta=zeros(1,nP);

mu=zeros(1,nP);
d=zeros(1,nP);
C=zeros(1,nP);
CG_ratio_d_constant=zeros(1,nP);
R0_d_C=zeros(1,nP);

mu_d_R=zeros(1,nP);
rho_d_R=zeros(1,nP);
CG_ratio_d_R=zeros(1,nP);
R0_d_R=zeros(1,nP);

mu_d_R2=zeros(1,nP);
rho_d_R2=zeros(1,nP);
CG_ratio_d_R2=zeros(1,nP);
R0_d_R2=zeros(1,nP);

mu_d_R3=zeros(1,nP);
rho_d_R3=zeros(1,nP);
CG_ratio_d_R3=zeros(1,nP);
R0_d_R3=zeros(1,nP);

mu_d_Rbeta=zeros(1,nP);
rho_d_Rbeta=zeros(1,nP);
CG_ratio_d_Rbeta=zeros(1,nP);
beta=zeros(1,nP);
R0_beta=zeros(1,nP);

a=zeros(1,nP);
b=zeros(1,nP);
V0_exp=zeros(1,nP);
%number of observations in each sample
nn=zeros(nP,1);

for i=1:nP
    pdv=data_v.(sprintf('v%d',i));
    pdt=data_time.(sprintf('t%d',i));

    n=length(pdv);
    nn(i,1)=n;


    %%%plolt patient data
    figure

    %do not fit initial values at beginning of experiment where
    %tumor volume is zero.
    plot(pdt(1:n), pdv(1:n), 'o')
    
    
    %%%%%Initial value for simulation
   
    r0=(pdv(1)*3/(4*pi))^(1/3);
    R0ub=r0;
    R0lb=r0;

   %suppose measurement error is almost always between -1 and 1 
   if strcmp(error_type,'meas_error')==1
       R0ub=r0+1;
       R0lb=r0-1;
       R0lb=max(R0lb,10^(-3));
   end
   
    
    %choose upper and lower bounds for parameters

    %%%set upper and lower bound for the estimate of mu
    mlb=0;
    %mlb=log(2)/1000;     
    mub=log(2);    
    %Set upper and lower bound for the estimate of d
    
    dlb=7*10^-3;
    dub=20;
    
   
    %Set upper and lower bound for the estimate of C
    Clb=10^5;
    Cub=3*10^6;
    %Set upper and lower bound for the estimate of CG_ratio
    CG_ratio_d_constantlb=1;
    CG_ratio_d_constantub=100;
    
    %select initial guess for parameters

    %initial guess for mu
    if (~exist('mu_d_constant0', 'var'))
        mu0=0.0019*6;
    else
        mu0=mu_d_constant0(i);
    end
    %initial guess for d
    if (~exist('d00', 'var'))
        d0=15*10^-3;
    else
        d0=d00(i);
    end
    %initial guess for C
    if (~exist('C00', 'var'))
        C0=10^5;
    else 
        C0=C00(i);
    end
    %initial guess for CG_ratio
    if (~exist('CG_ratio_d_constant00', 'var'))
        CG_ratio_d_constant0=10;
    else
        CG_ratio_d_constant0=CG_ratio_d_constant00(i);
    end
    %initial guess for R0 if error is due to measurement
    if strcmp(error_type,'meas_error')==1
    if (~exist('R0_d_C00', 'var'))
        R0_d_C0=r0;
    else
        R0_d_C0=R0_d_C00(i);
    end
    else
        R0_d_C0=r0;
    end
    
    %Set upper and lower bound vector 
    lb=[mlb dlb Clb CG_ratio_d_constantlb,R0lb];
    ub=[mub dub Cub  CG_ratio_d_constantub,R0ub];
    
    %initial guess vector 
    [pp0]=[mu0 d0 C0  CG_ratio_d_constant0,R0_d_C0];
    
    %%%%estimate the parameters
    SE=@(pp)SE_fun(pp,pdt,pdv,error_type);
    options = optimoptions(@fmincon,'MaxFunctionEvaluations',100000,'MaxIterations',10000);
    [pp,error_d_constant(i)] = fmincon(SE, pp0,[],[],[],[],lb,ub,[],options);
    mu(i)=pp(1);
    d(i)=pp(2);
    C(i)=pp(3);
    CG_ratio_d_constant(i)=pp(4);
    R0_d_C(i)=pp(5);
    
    fprintf('model d=C, sample=%0.0f\n',i)
    fprintf('Estimated mu is %9.4f\n', mu(i))
    fprintf('Estimated d is %9.4f\n', d(i))
    fprintf('Estimated C is %9.4f\n', C(i))
    fprintf('Estimated C:G ratio is %9.4f\n',CG_ratio_d_constant(i))
    
    %%%%use the estimated paprameters to graph the solutions 
    
    tumor_growth_model_fun=@(t,y)tumor_growth_model_dC(t,y,mu(i),d(i),C(i),CG_ratio_d_constant(i));    %%%%%%%%
    [t,y] = ode45(tumor_growth_model_fun,[pdt(1) max(pdt)],R0_d_C(i));
    fitr=y(:,1);
    fitv=(4*pi*fitr.^3)/3;
    %plot solution time series
    hold on
    plot(t, fitv, 'LineWidth', 1,'Color','r')
    ylabel('volume (mm$^3$)','Interpreter','latex','FontSize',12)
    xlabel('days','Interpreter','latex','FontSize',12)
    
%     %compute and store residuals
%     [~,y] = ode45(tumor_growth_model_fun,pdt(1:n),R0_d_C(i));
%     fitr=y(:,1);
%     fitv=(4*pi*fitr.^3)/3;
%     
%     res_d_constant.(sprintf('v%d',i))=pdv(1:n)-fitv;

    %plot solution increments
    pdr=((3/(4*pi))*pdv).^(1/3);
    R=zeros(n,1);
    for ii=2:n
        [~,y]=ode45(tumor_growth_model_fun,[pdt(ii-1) pdt(ii)],pdr(ii-1));
        R(ii)=y(length(y));
    end
    fitV=(4*pi*R.^3)/3;
    plot(pdt(2:n),fitV(2:n),'*','Color','r');


    %%%% Fit the allometric growth model d = rho*R%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %Set upper and lower bound for the estimate of rho
    rholb=0;
    rhoub=10;
    %Set upper and lower bound for the estimate of mu
    %mlb=log(2)/100;
    mlb=0;
    %mub=log(2);
    mub=log(2);
    %Set upper and lower bound for the estimate of CG_ratio.
   
    CG_ratiolb=20;
    CG_ratioub=inf;

    %select initial guess for parameters
    %guess for rho
    
    rho0=0.1;
    
    %guess for mu
    mu0=0.0019*6;
    
    %guess for CG_ratio
    CG_ratio0=40;
    
    %initial guess for R0 if error is due to measurement
    R0=r0;
    
    %initial vector 
    [ppallo0]=[mu0 rho0 CG_ratio0, R0];
    
    %Set upper and lower bound vector 
    lb=[mlb rholb CG_ratiolb R0lb];
    ub=[mub rhoub CG_ratioub R0ub];
    
    SE=@(ppallo)SE_fun_allometric_d_R(ppallo,pdt,pdv,error_type);
    [ppallo,error_d_R(i)] = fmincon(SE,ppallo0,[],[],[],[],lb,ub,[],options);
    mu_d_R(i)=ppallo(1);
    rho_d_R(i)=ppallo(2);
    CG_ratio_d_R(i)=ppallo(3);
    R0_d_R(i)=ppallo(4);
    
    fprintf('model d = rho R sample=%0.0f\n',i)
    fprintf('Estimated mu is %9.4f\n', mu_d_R(i))
    fprintf('Estimated rho is %9.4f\n', rho_d_R(i))
    fprintf('Estimated CG_ratio is %9.4f\n', CG_ratio_d_R(i))
    
    %%%%use the estimated parameters to graph the solutions 
    
    
%     fun=@(t,y)tumor_allometric_growth_model_d_Rbeta(t,y,mu_d_R(i),rho_d_R(i),CG_ratio_d_R(i),1);    %%%%%%%%
%     [t,y] = ode45(fun,[pdt(1) max(pdt)],R0_d_R(i));
%     %plot solution as a time series
%     fitr=y(:,1);
%     fitv=(4*pi*fitr.^3)/3;
%     plot(t, fitv, 'LineWidth', 1)
% 
%     %compute and store residuals
%     [~,y] = ode45(fun,pdt(1:n),R0_d_R(i));
%     fitr=y(:,1);
%     fitv=(4*pi*fitr.^3)/3;
%     %plot solution increments
%     pdr=((3/(4*pi))*pdv).^(1/3);
%     R=zeros(n,1);
%     for ii=2:n
%         [~,y]=ode45(fun,[pdt(ii-1) pdt(ii)],pdr(ii-1));
%         R(ii)=y(length(y));
%     end
%     fitV=(4*pi*R.^3)/3;
%     plot(pdt(2:n),fitV(2:n),'*');
    
    
%     res_d_R.(sprintf('v%d',i))=pdv(1:n)-fitv';
    
    %%%% Fit the allometric growth model d=rho*R^2%%%%%%%%%%%%%%%%%%%%%%%%

    SE=@(ppallo)SE_fun_allometric_d_R2(ppallo,pdt,pdv,error_type);
    [ppallo,error_d_R2(i)] = fmincon(SE,ppallo0,[],[],[],[],lb,ub,[],options);
    mu_d_R2(i)=ppallo(1);
    rho_d_R2(i)=ppallo(2);
    CG_ratio_d_R2(i)=ppallo(3);
    R0_d_R2(i)=ppallo(4);
    
    fprintf('model d = rho R^2, sample=%0.0f\n',i)
    fprintf('Estimated mu is %9.4f\n', mu_d_R2(i))
    fprintf('Estimated rho is %9.4f\n', rho_d_R2(i))
    fprintf('Estimated CG_ratio is %9.4f\n', CG_ratio_d_R2(i))
    
    %%%%use the estimated paprameters to graph the solutions 
    
    
    
%     fun=@(t,y)tumor_allometric_growth_model_d_Rbeta(t,y,mu_d_R2(i),rho_d_R2(i),CG_ratio_d_R2(i),2);    %%%%%%%%
%     [t,y] = ode45(fun,[pdt(1) max(pdt)],R0_d_R2(i));
%     fitr=y(:,1);
%     fitv=(4*pi*fitr.^3)/3;
%     hold on
%     %plot solution as a time series
%     plot(t, fitv, 'LineWidth', 1)
% 
% 
%     %compute and store residuals
%     [~,y] = ode45(fun,pdt(1:n),R0_d_R2(i));
%     fitr=y(:,1);
%     fitv=(4*pi*fitr.^3)/3;
%     %compute residuals
%     res_d_R2.(sprintf('v%d',i))=pdv(1:n)-fitv';
    
    %%%% Fit the allometric growth model d = rho*R^3  %%%%%%%%
    %%%%%%%%% with undetected low density growth region %%%%%%
    
    
    SE=@(ppallo)SE_fun_allometric_d_R3(ppallo,pdt,pdv,error_type);
    [ppallo,error_d_R3(i)] = fmincon(SE,ppallo0,[],[],[],[],lb,ub,[], options);
    mu_d_R3(i)=ppallo(1);
    rho_d_R3(i)=ppallo(2);
    CG_ratio_d_R3(i)=ppallo(3);
    R0_d_R3(i)=ppallo(4);
    
    fprintf('model d = rho R^3,sample=%0.0f\n',i)
    fprintf('Estimated mu is %9.4f\n', mu_d_R3(i))
    fprintf('Estimated rho is %9.4f\n', rho_d_R3(i))
    fprintf('Estimated CG_ratio is %9.4f\n', CG_ratio_d_R3(i))
    
    %%%%use the estimated paprameters to graph the solutions 
    
    
    
%     fun=@(t,y)tumor_allometric_growth_model_d_Rbeta(t,y,mu_d_R3(i),rho_d_R3(i),CG_ratio_d_R3(i),3);    %%%%%%%%
%     [t,y] = ode45(fun,[pdt(1) max(pdt)],R0_d_R3(i));
%     fitr=y(:,1);
%     fitv=(4*pi*fitr.^3)/3;
%     %plot solution as a time series
%     hold on
%     plot(t, fitv, 'LineWidth', 1)
%     
%     
%     %compute and store residuals
%     [~,y] = ode45(fun,pdt(1:n),R0_d_R3(i));
%     fitr=y(:,1);
%     fitv=(4*pi*fitr.^3)/3;
%     
%     res_d_R3.(sprintf('v%d',i))=pdv(1:n)-fitv;


     %%%% Fit the allometric growth model d = rho*R^beta  %%%%%%%%
    %%%%%%%%% with undetected low density growth region %%%%%%

    %set upper and lower bounds for estimate of beta
    betalb=0;
    betaub=11;
    %betaub=3;
    %select initial guess for parameters
    %guess for rho
    if (~exist('rho_d_Rbeta0', 'var'))
        rho0=0.1;
    else
        rho0=rho_d_Rbeta0(i);
    end
    %guess for mu
    if (~exist('mu_d_Rbeta0', 'var'))
        mu0=0.0019*6;
    else
        mu0=mu_d_Rbeta0(i);
    end
    %guess for CG_ratio
    if (~exist('CG_ratio_d_Rbeta0', 'var'))
        CG_ratio0=40;
    else
        CG_ratio0=CG_ratio_d_Rbeta0(i);
    end
    %guess for beta
    if (~exist('beta00', 'var'))
        beta0=1.2;
    else
        beta0=beta00(i);
    end

     %initial guess for R0 if error is due to measurement
    if strcmp(error_type,'meas_error')==1
    if (~exist('R0_beta00', 'var'))
        R0_beta0=r0;
    else
        R0_beta0=R0_beta00(i);
    end
    else
        R0_beta0=r0;
    end
    
    
    %initial vector 
    [ppallo0]=[mu0 rho0 CG_ratio0, beta0,R0_beta0];
    
    %Set upper and lower bound vector 
    lb=[mlb rholb CG_ratiolb betalb R0lb];
    ub=[mub rhoub CG_ratioub betaub R0ub];
    
    SE=@(ppallo)SE_fun_allometric_d_beta(ppallo,pdt,pdv,error_type);
    [ppallo,error_d_Rbeta(i)] = fmincon(SE,ppallo0,[],[],[],[],lb,ub,[],options);
    
    %see if the fit is better than the other models in the same class
    [min_lse, min_index]=min([error_d_R(i) error_d_R2(i) error_d_R3(i)]);
    if min_lse<error_d_Rbeta(i)
        if min_index==1
        [ppallo0]=[mu_d_R(i) rho_d_R(i) CG_ratio_d_R(i), 1, R0_d_R(i)];
        end
        if min_index==2
        [ppallo0]=[mu_d_R2(i) rho_d_R2(i) CG_ratio_d_R2(i), 2, R0_d_R2(i)];
        end
        if min_index==3
        [ppallo0]=[mu_d_R3(i) rho_d_R3(i) CG_ratio_d_R3(i), 3, R0_d_R3(i)];
        end

    SE=@(ppallo)SE_fun_allometric_d_beta(ppallo,pdt,pdv,error_type);
    [ppallo_new,lse_d_Rbeta_new] = fmincon(SE,ppallo0,[],[],[],[],lb,ub,[],options);

    if lse_d_Rbeta_new<error_d_Rbeta(i)
        ppallo=ppallo_new;
        error_d_Rbeta(i)=lse_d_Rbeta_new;
    end
    end

    mu_d_Rbeta(i)=ppallo(1);
    rho_d_Rbeta(i)=ppallo(2);
    CG_ratio_d_Rbeta(i)=ppallo(3);
    beta(i)=ppallo(4);
    R0_beta(i)=ppallo(5);
    
    fprintf('model d = rho R^beta, sample=%0.0f\n',i)
    fprintf('Estimated mu is %9.4f\n', mu_d_Rbeta(i))
    fprintf('Estimated rho is %9.4f\n', rho_d_Rbeta(i))
    fprintf('Estimated CG_ratio is %9.4f\n', CG_ratio_d_Rbeta(i))
    fprintf('Estimated beta is %9.4f\n', beta(i))
    
    %%%%use the estimated paprameters to graph the solutions 
    
    fun=@(t,y)tumor_allometric_growth_model_d_Rbeta(t,y,mu_d_Rbeta(i),rho_d_Rbeta(i),CG_ratio_d_Rbeta(i),beta(i));    %%%%%%%%
    [t,y] = ode45(fun,[pdt(1) max(pdt)],R0_beta(i));
    fitr=y(:,1);
    fitv=(4*pi*fitr.^3)/3;
    %plot solution as a time series
    hold on
    plot(t, fitv, 'LineWidth', 1,'Color','g')
    
    
    %compute and store residuals
    [~,y] = ode45(fun,pdt(1:n), R0_beta(i));
    fitr=y(:,1);
    fitv=(4*pi*fitr.^3)/3;
    
    res_d_beta.(sprintf('v%d',i))=pdv(1:n)-fitv;
    
    %plot solution increments
    pdr=((3/(4*pi))*pdv).^(1/3);
    R=zeros(n,1);
    for ii=2:n
        [~,y]=ode45(fun,[pdt(ii-1) pdt(ii)],pdr(ii-1));
        R(ii)=y(length(y));
    end
    fitV=(4*pi*R.^3)/3;
    plot(pdt(2:n),fitV(2:n),'*','Color','g');
    
    
    
    %%% Fit Exponential model %%%%%%%

    %set upper and lower bounds for b
    blb=0;
    bub=2;

    %set upper and lower bounds for V0
    %suppose measurement error is almost always between -1 and 1 
   if strcmp(error_type,'meas_error')==1
       V0ub=(4*pi*R0ub^3)/3;
       V0lb=(4*pi*R0lb^3)/3;
   else
       V0ub=pdv(1);
       V0lb=pdv(1);
   end
    %initial guess for b
    if (~exist('b00', 'var'))
        b0=1.2;
    else
        b0=b00(i);
    end

     %initial guess for V0 if error is due to measurement
    if strcmp(error_type,'meas_error')==1
    if (~exist('V0_exp00', 'var'))
        V0_exp0=pdv(1);
    else
        V0_exp0=V0_exp00(i);
    end
    else
        V0_exp0=pdv(1);
    end
    
    
    %find the optimal value of b
    SE=@(b)SE_fun_exp_modelb(b,pdt,pdv,error_type);
    [x,error_exponential(i)] = fmincon(SE, [b0 V0_exp0], [], [], [], [],[blb V0lb],[bub V0ub],[],options);
    b(i)=x(1);
    V0_exp(i)=x(2);
    
    %compute the cooresponding optimal value of a
    
    if b(i)>1
        aub=(V0_exp(i)^(1-b(i)))/((b(i)-1)*max(pdt));
        %stipulate that per capita division rate in largest observed tumor
        %is not greater than 1 division per day
        aub=min(aub,(max(pdv)^(1-b(i)))*log(2));
        alb=0;
        a0=(alb+aub)/2;
    end
    
    if b(i)<=1
        %aub=[];
        %stipulate per capita division rate in nacent tumor is at most one
        %division per day
        aub=log(2)*10^(-5);
        alb=0;
        a0=1;
    end 
    SE=@(a)SE_fun_exp_model([a b(i) V0_exp(i)],pdt,pdv,error_type);
    [a(i),~] = fmincon(SE, a0, [], [], [], [],alb,aub,[],options);
    aa=a(i);
    bb=b(i);
    
    fprintf('model dV/dt = alpha V^beta,sample=%0.0f\n',i)
    fprintf('Estimated a is %9.4f\n', aa)
    fprintf('Estimated b is %9.4f\n', bb)
    
    %%%%use the estimated paprameters to graph the solutions 
    
    if bb==1
        simv_exp=@(t)V0_exp(i)*exp(aa*(t-pdt(1)));
    else
        simv_exp= @(t)(aa*(1-bb)*(t-pdt(1)) + V0_exp(i)^(1-bb)).^(1/(1-bb));
    end
    
    v_exp=simv_exp(t);
    %plot solution time series
    hold on
    plot(t, v_exp, 'LineWidth', 1,'Color','b')
    

    %plot solution increments
    if bb==1
            simv_exp=@(t,i)pdv(i)*exp(aa*(t-pdt(i)));
        else
            simv_exp= @(t,i)(aa*(1-bb)*(t-pdt(i)) + pdv(i)^(1-bb)).^(1/(1-bb));
    end

    for ii=2:n
        fitV(ii)=simv_exp(pdt(ii),ii-1);
    end
    plot(pdt(2:n),fitV(2:n),'*','Color','b');
    
    %compute and store residuals
    
%     fitv=simv_exp(pdt(1:n));
%     res_exp.(sprintf('v%d',i))=pdv(1:n)-fitv;
    
    %annotate plot
    
    %add legend
    
    
    lbl1= '$d=C$';
    lbl2= '$d=C$';
    lbl3= '$d=\\rho R^{\\beta}$';
    lbl4= '$d=\\rho R^{\\beta}$';
    lbl5= '$dV/dt=aV^{b}$';
    lbl6= '$dV/dt=aV^{b}$';
   
    
    legend('patient data',sprintf(lbl1),sprintf(lbl2),sprintf(lbl3),sprintf(lbl4),sprintf(lbl5),sprintf(lbl6),'Interpreter','latex','FontSize',10,'Location','northwest')
    
    
    ttl=sprintf(' (sample %d)',i);
    ttl=strcat(data_source,ttl);
    title(ttl,'Interpreter','latex');

    figname=sprintf('_%d.pdf',i);
    figname=strcat(data_source,figname);
    figname=strcat(error_type,figname);
    exportgraphics(gcf,figname);

    hold off
end

%plot residuals for exponential model
% for i=1:nP
%     hold off
%     figure
%     t=data_time.(sprintf('t%d',i));
%     res=res_exp.(sprintf('v%d',i));
%     plot(t, res, 'o');
%     ylabel('residual ($\mathrm{mm}^3$)','Interpreter','latex','FontSize',12)
%     xlabel('days','Interpreter','latex','FontSize',12)
%     title(sprintf('exponential model residuals, sample %d',i),'Interpreter','latex')
%     exportgraphics(gcf,sprintf('exp_res_%d.eps',i));
% end

%create matrix of lse for ranking

error_matrix=[error_d_constant' error_d_Rbeta' error_exponential'];
[~,sample_rank_index]=sort(error_matrix, 2);
sample_rank=zeros(nP,3);
%get ranking of each model
for i=1:nP
    for j=1:3
sample_rank(i,j) = find(sample_rank_index(i,:)==j);
    end
end

% lse_matrix=[lse_d_constant' lse_d_R' lse_d_R2' lse_d_R3' lse_d_Rbeta' lse_exponential'];
% [~,sample_rank_index]=sort(lse_matrix, 2);
% sample_rank=zeros(nP,6);
% %get ranking of each model
% for i=1:nP
%     for j=1:6
% sample_rank(i,j) = find(sample_rank_index(i,:)==j);
%     end
% end

%create, print, and export tables of errors and parameter estimates

filename = 'parameter_fits.xlsx';
filename=strcat(data_source,filename);
filename=strcat(error_type,filename);

table_d_constant_bf_parameters=table(mu',d',C',CG_ratio_d_constant',sample_rank(:,1),'VariableNames', {'mu','d','C','CG ratio','rank'})
writetable(table_d_constant_bf_parameters,filename,'Sheet',1)


figurename='bf_parameters_d_constant.pdf';
figurename=strcat(data_source,figurename);
figurename=strcat(error_type,figurename);
fig = uifigure;
uitable(fig,'Data',table_d_constant_bf_parameters{:,:},'ColumnWidth','fit','ColumnName',table_d_constant_bf_parameters.Properties.VariableNames,'Units', 'Normalized', 'Position',[0 0 1 1])
pause(5);

exportapp(fig,figurename);
% 
% 
% table_d_R_bf_parameters=table(mu_d_R',rho_d_R',CG_ratio_d_R',sample_rank(:,2),'VariableNames', {'mu_d_R','rho_d_R','CG_ratio_d_R','sample_rank'})
% writetable(table_d_R_bf_parameters,filename,'Sheet',2)
% 
% figurename='bf_parameters_d_R.pdf';
% figurename=strcat(data_source,figurename);
% fig = uifigure;
% uitable(fig,'Data',table_d_R_bf_parameters{:,:},'ColumnWidth','fit','ColumnName',table_d_R_bf_parameters.Properties.VariableNames,'Units', 'Normalized', 'Position',[0 0 1 1])
% exportapp(fig,figurename);
% 
% 
% table_d_R2_bf_parameters=table(mu_d_R2',rho_d_R2',CG_ratio_d_R2',sample_rank(:,3),'VariableNames', {'mu_d_R2','rho_d_R2','CG_ratio_d_R2','sample_rank'})
% writetable(table_d_R2_bf_parameters,filename,'Sheet',3)
% 
% figurename='bf_parameters_d_R2.pdf';
% figurename=strcat(data_source,figurename);
% fig = uifigure;
% uitable(fig,'Data',table_d_R2_bf_parameters{:,:},'ColumnWidth','fit','ColumnName',table_d_R2_bf_parameters.Properties.VariableNames,'Units', 'Normalized', 'Position',[0 0 1 1])
% exportapp(fig,figurename);
% 
% table_d_R3_bf_parameters=table(mu_d_R3',rho_d_R3',CG_ratio_d_R3',sample_rank(:,4),'VariableNames', {'mu_d_R3','rho_d_R3','CG_ratio_d_R3','sample_rank'})
% writetable(table_d_R3_bf_parameters,filename,'Sheet',4)
% 
% figurename='bf_parameters_d_R3.pdf';
% figurename=strcat(data_source,figurename);
% fig = uifigure;
% uitable(fig,'Data',table_d_R3_bf_parameters{:,:},'ColumnWidth','fit','ColumnName',table_d_R3_bf_parameters.Properties.VariableNames,'Units', 'Normalized', 'Position',[0 0 1 1])
% exportapp(fig,figurename);

table_d_Rbeta_bf_parameters=table(mu_d_Rbeta',rho_d_Rbeta',CG_ratio_d_Rbeta', beta',sample_rank(:,2),'VariableNames', {'mu','rho','CG ratio','beta','rank'})
writetable(table_d_Rbeta_bf_parameters,filename,'Sheet',2)

figurename='bf_parameters_d_Rbeta.pdf';
figurename=strcat(data_source,figurename);
figurename=strcat(error_type,figurename);
fig = uifigure;
uitable(fig,'Data',table_d_Rbeta_bf_parameters{:,:},'ColumnWidth','fit','ColumnName',table_d_Rbeta_bf_parameters.Properties.VariableNames,'Units', 'Normalized', 'Position',[0 0 1 1])
pause(5);
exportapp(fig,figurename);

table_exponential_bf_parameters=table(a',b',sample_rank(:,3),'VariableNames', {'a','b','rank'})
writetable(table_exponential_bf_parameters,filename,'Sheet',3)

figurename='bf_parameters_exp.pdf';
figurename=strcat(data_source,figurename);
figurename=strcat(error_type,figurename);
fig = uifigure;
uitable(fig,'Data',table_exponential_bf_parameters{:,:},'ColumnWidth','fit','ColumnName',table_exponential_bf_parameters.Properties.VariableNames,'Units', 'Normalized', 'Position',[0 0 1 1])
pause(5);
exportapp(fig,figurename);

table_lse=table(error_d_constant', error_d_Rbeta', error_exponential','VariableNames', {'model 1 (d=C)', 'model 2 (d=R^beta)', 'model 3 (dV/dt=aV^b)'})
writetable(table_lse,filename,'Sheet',7)

figurename='lse.pdf';
figurename=strcat(data_source,figurename);
figurename=strcat(error_type,figurename);
fig = uifigure;
uit=uitable(fig,'Data',table_lse{:,:},'ColumnWidth','fit','ColumnName',table_lse.Properties.VariableNames,'Units', 'Normalized', 'Position',[0 0 1 1])

[~,minIndices] = min(error_matrix,[],2);
minIndices=[(1:nP)', minIndices];

s = uistyle('BackgroundColor',[1 0.6 0.6]);
addStyle(uit,s,'cell',minIndices);

pause(5);
exportapp(fig,figurename);



if strcmp(error_type,'meas_error')==1
    [~,minIndices_relL] = min(error_matrix(:,1:2),[],2);
    relL=zeros(nP,2);
    num_param=zeros(1,2);
    num_param(1)=5;
    num_param(2)=5;
    for k=1:nP
        for j=1:2
            relL(k,j)=exp(num_param(minIndices_relL(k))-num_param(j)+(error_matrix(k,minIndices_relL(k))-error_matrix(k,j))/(2*.01));
        end
    end
table_relL=table(relL(:,1), relL(:,2),'VariableNames', {'model 1 (d=C)', 'model 2 (d=R^beta)'})

figurename='relL.pdf';
figurename=strcat(data_source,figurename);
figurename=strcat(error_type,figurename);
fig = uifigure;
uit=uitable(fig,'Data',table_relL{:,:},'ColumnWidth','fit','ColumnName',table_relL.Properties.VariableNames,'Units', 'Normalized', 'Position',[0 0 1 1])
pause(5);
exportapp(fig,figurename);
end

if strcmp(error_type,'meas_error')==1
    AICc=zeros(nP,3);
    num_param=zeros(1,3);
    num_param(1)=5;
    num_param(2)=5;
    num_param(3)=3;
    for k=1:nP
        for j=1:3
            if nn(k)<=num_param(j)
                AICc(k,j)=inf;
            else
                AICc(k,j)=2*(num_param(j)+(error_matrix(k,j)/(2*.01)))+2*((num_param(j))^2+num_param(j))/(nn(k)-num_param(j)-1);
            end
        end
    end
    [~,minIndices_relL3] = min(AICc,[],2);
    relL3=zeros(nP,3);
    for k=1:nP
        for j=1:3
            relL3(k,j)=exp((AICc(k,minIndices_relL3(k))-AICc(k,j))/2);
        end
    end
    table_AICc=table(AICc(:,1),AICc(:,2),AICc(:,3),'VariableNames', {'model 1 (d=C)', 'model 2 (d=R^beta)','model 3 (dV/dt=aV^b)'})

    table_relL3=table(relL3(:,1),relL3(:,2),relL3(:,3),'VariableNames', {'model 1 (d=C)', 'model 2 (d=R^beta)','model 3 (dV/dt=aV^b)'})

    figurename='relL3.pdf';
    figurename=strcat(data_source,figurename);
    figurename=strcat(error_type,figurename);
    fig = uifigure;
    uit=uitable(fig,'Data',table_relL3{:,:},'ColumnWidth','fit','ColumnName',table_relL3.Properties.VariableNames,'Units', 'Normalized', 'Position',[0 0 1 1])
    pause(5);
    exportapp(fig,figurename);
end



end

%functions for computing error between model simulations and data

function error=SE_fun(x, pdt, pdv,error_type)  
%%%for a particular parameter set 'x', this function computes the error
%between the simulated data for the growth region model d=C and real tumor volume data

%%%%%The following are parameters to be estimated
mu=x(1);
d=x(2);
C=x(3);
CG_ratio=x(4);
R0=x(5);

%%%%Solve the model and compute the error
fun=@(t,y)tumor_growth_model_dC(t,y,mu,d,C,CG_ratio);%%%%%%%%
n=length(pdt);
simv=zeros(n,1);
simr=zeros(n,1);
if strcmp(error_type, 'LL_error')~=1
[~,y] = ode45(fun,pdt(1:n),R0);
simr=y(:,1);
simv=(4*pi*simr.^3)/3;
end

if strcmp(error_type, 'LL_error')==1
    pdr=((3/(4*pi))*pdv).^(1/3);
    for i=2:n
        [~,y]=ode45(fun,[pdt(i-1) pdt(i)],pdr(i-1));
        simr(i)=y(length(y));
    end
   simv=(4*pi*simr.^3)/3;
end

%return absolute or relative error
[error]=error_fun(pdt,pdv,simr,simv,error_type);      
end

function error=SE_fun_allometric_d_R(x,pdt, pdv,error_type)  
%%%for a particular parameter set 'x', this function computes the error
%between the model d=R and real tumor volume data

%%%%%The following are parameters to be estimated
mu=x(1);
rho=x(2);
CG_ratio=x(3);
R0=x(4);

%%%%Solve the model and compute the error
fun=@(t,y)tumor_allometric_growth_model_d_Rbeta(t,y,mu,rho,CG_ratio,1);    %%%%%%%%
n=length(pdt);
simv=zeros(n,1);
simr=zeros(n,1);
if strcmp(error_type, 'LL_error')~=1
[~,y] = ode45(fun,pdt(1:n),R0);
simr=y(:,1);
simv=(4*pi*simr.^3)/3;
end

if strcmp(error_type, 'LL_error')==1
    pdr=((3/(4*pi))*pdv).^(1/3);
    for i=2:n
        [~,y]=ode45(fun,[pdt(i-1) pdt(i)],pdr(i-1));
        simr(i)=y(length(y));
    end
   simv=(4*pi*simr.^3)/3;
end
%return absolute or relative error
[error]=error_fun(pdt,pdv,simr,simv,error_type);  
end


function error=SE_fun_allometric_d_R2(x,pdt, pdv,error_type)  
%%%for a particular parameter set 'x', this function computes the error
%between the mode d=R2 and real tumor volume data

%%%%%The following are parameters to be estimated
mu=x(1);
rho=x(2);
CG_ratio=x(3);
R0=x(4);

%%%%Solve the model and compute the error
fun=@(t,y)tumor_allometric_growth_model_d_Rbeta(t,y,mu,rho,CG_ratio,2);    %%%%%%%%
n=length(pdt);
simv=zeros(n,1);
simr=zeros(n,1);
if strcmp(error_type, 'LL_error')~=1
[~,y] = ode45(fun,pdt(1:n),R0);
simr=y(:,1);
simv=(4*pi*simr.^3)/3;
end

if strcmp(error_type, 'LL_error')==1
    pdr=((3/(4*pi))*pdv).^(1/3);
    for i=2:n
        [~,y]=ode45(fun,[pdt(i-1) pdt(i)],pdr(i-1));
        simr(i)=y(length(y));
    end
   simv=(4*pi*simr.^3)/3;
end

%return absolute or relative error
[error]=error_fun(pdt,pdv,simr,simv,error_type);  
end

function error=SE_fun_allometric_d_R3(x, pdt, pdv,error_type)  
%%%for a particular parameter set 'x', this function computes the error
%between the mode d=R^3 and real tumor volume data

%%%%%The following are parameters to be estimated
mu=x(1);
rho=x(2);
CG_ratio=x(3);
R0=x(4);

%%%%Solve the model and compute the error
fun=@(t,y)tumor_allometric_growth_model_d_Rbeta(t,y,mu,rho,CG_ratio,3);    %%%%%%%%
n=length(pdt);
simv=zeros(n,1);
simr=zeros(n,1);
if strcmp(error_type, 'LL_error')~=1
[~,y] = ode45(fun,pdt(1:n),R0);
simr=y(:,1);
simv=(4*pi*simr.^3)/3;
end

if strcmp(error_type, 'LL_error')==1
    pdr=((3/(4*pi))*pdv).^(1/3);
    for i=2:n
        [~,y]=ode45(fun,[pdt(i-1) pdt(i)],pdr(i-1));
        simr(i)=y(length(y));
    end
   simv=(4*pi*simr.^3)/3;
end

%return absolute or relative error
[error]=error_fun(pdt,pdv,simr,simv,error_type);
end

function error=SE_fun_allometric_d_beta(x,pdt,pdv,error_type)
%%%for a particular parameter set 'x', this function computes the error
%between the model d=R^beta and real tumor volume data

%%%%%The following are parameters to be estimated
mu=x(1);
rho=x(2);
CG_ratio=x(3);
beta=x(4);
R0=x(5);

%%%%Solve the model and compute the error
fun=@(t,y)tumor_allometric_growth_model_d_Rbeta(t,y,mu,rho,CG_ratio,beta);    %%%%%%%%
nn=length(pdt);
simv=zeros(nn,1);
simr=zeros(nn,1);
if strcmp(error_type, 'LL_error')~=1
[~,y] = ode45(fun,pdt(1:nn),R0);
simr=y(:,1);
simv=(4*pi*simr.^3)/3;
end

if strcmp(error_type, 'LL_error')==1
    pdr=((3/(4*pi))*pdv).^(1/3);
    for i=2:nn
        [~,y]=ode45(fun,[pdt(i-1) pdt(i)],pdr(i-1));
        simr(i)=y(length(y));
    end
   simv=(4*pi*simr.^3)/3;
end


%return absolute or relative error
[error]=error_fun(pdt,pdv,simr,simv,error_type);    
end



function error=SE_fun_exp_model(x,pdt,pdv,error_type) 
%%%for a particular parameter set 'x', this function computes the error
%between the simulation and real tumor volume data

%%%%%The following are parameters to be estimated
a=x(1);
b=x(2);
V0=x(3);
n=length(pdt);
simv=zeros(n,1);

%%%%Solve the model and compute the error
if strcmp(error_type, 'LL_error')~=1
    if b==1
        simv=V0*exp(a*(pdt(1:n)-pdt(1)));
    end
    if b~=1
        simv= (a*(1-b)*(pdt(1:n)-pdt(1)) + V0^(1-b)).^((1)/(1-b));
    end
end

if strcmp(error_type, 'LL_error')==1
    if b==1
        for ii=2:n
            simv(ii)=pdv(ii-1)*exp(a*(pdt(ii)-pdt(ii-1)));
        end
    end
    if b~=1
        for ii=2:n
            simv(ii)= (a*(1-b)*(pdt(ii)-pdt(ii-1)) + pdv(ii-1)^(1-b)).^((1)/(1-b));
        end
    end
end

simr=(3*simv/(4*pi)).^(1/3);

%return absolute or relative error
[error]=error_fun(pdt,pdv,simr,simv,error_type);
if error==inf
    error=realmax;
end
end



function lse_exponential=SE_fun_exp_modelb(x,pdt,pdv,error_type) 
%%%for a particular parameter set 'x', this function computes the error
%between the simulation and real tumor volume data
%compute upper and lower bounds for a
b=x(1);
V0=x(2);
nn=length(pdt);
if b>1
    if strcmp(error_type,'LL_error')~=1
        aub=(V0^(1-b))/((b-1)*(max(pdt)-pdt(1)));
        %stipulate that per capita division rate in largest observed tumor
        %is not greater than 1 division per day
        aub=min(aub,(max(pdv)^(1-b))*log(2));
    end
    if strcmp(error_type,'LL_error')==1
        aub=inf;
        for ii=2:nn
            aub_new=(pdv(ii-1)^(1-b))/((b-1)*((pdt(ii))-pdt(ii-1)));
            aub=min(aub,aub_new);
            %stipulate that per capita division rate in largest observed tumor
            %is not greater than 1 division per day
            aub=min(aub,(max(pdv)^(1-b))*log(2));
        end
    end
alb=0;
a0=(alb+aub)/2;
end

if b<=1
    %aub=[];
    %stipulate per capita division rate in nacent tumor is at most one
    %division per day
    aub=log(2)*10^(-5);
    %stipulate that per capita division rate in unit tumor is not
    %greater than 1 division a day
    %aub=log(2);
    alb=0;
    a0=1;
end 

SE=@(a)SE_fun_exp_model([a b V0],pdt,pdv,error_type);
options = optimoptions(@fmincon,'MaxFunctionEvaluations',100000,'MaxIterations',10000,'Display','none');
[a,lse_exponential] = fmincon(SE, a0, [], [], [], [],alb,aub,[],options);

end

%function for measuring error
function [error]=error_fun(pdt,pdv,simr,simv,error_type)
if strcmp(error_type,'abs_error')==1
    error=sum((pdv-simv).^2);
end

if strcmp('rel_error', error_type)==1
    error=sum(((pdv-simv)./simv).^2);
end

if strcmp(error_type,'meas_error')==1
    pdr=(3*pdv/(4*pi)).^(1/3);
    error=sum((pdr-simr).^2);
end

if strcmp(error_type,'LL_error')==1
    n=length(pdt);
    tau=zeros(n,1);
    for i=2:n
    tau(i)=pdt(i)-pdt(i-1);
    end
    error=sum(abs(pdv(2:n)-simv(2:n))./(simv(2:n).*sqrt(tau(2:n))));
end

end