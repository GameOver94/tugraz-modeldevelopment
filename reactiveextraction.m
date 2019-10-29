% Model development
% WS 2019
%
%
%
% Reaktive Extration


clear variables; 
clc;
close all;
%%

% Liquid-liquid equilibrium
% x_i_p*gamma_i_p=x_i_dp*gamma_i_dp for i =1,2
% x_1 = 1- sum(x_i)
%n_1_p/n_1_dp=n_2_p/n_2_dp

alpha = 0.2;

options = optimoptions('fsolve','StepTolerance',1e-10);

%% calculation LLE

% critical Temperatur and composition

c=fsolve(@(x)critVal(x,alpha),[0.1,800]);

T_crit = c(2);
x_crit = c(1);

%% LLE diagram
i=1;
T(i)=T_crit;
for i=1:40
           
    
        %h=fsolve(@(x)dge(x,T(i),alpha),[x_crit-.2;x_crit+.2]);
        %min_cg(fun, x0, gtol, eps, maxiter, rho, delta, mu)
        h=min_cg(@(x)target(x,T(i),alpha),[x_crit-.2;x_crit+.2], 10^-7, 0.0001, 100000, 0.6, 0.01, 0.1);
        x1(i) = h(1);
        x1p(i) = h(2);
        
        temp=T;
         T(i+1) = T(i)-5;
         
end
    
   plot(x_crit,T_crit,'b--o',x1,temp,'b--o',x1p,temp,'b--o')
   xlim([0 1]);
   title ('Water-hexane');
   ylabel('temperature (K)');
   xlabel('mol fraction');

   
   %% tau calculation
   function tij = tau(T,i)

if i==1 %parameter ij
    A=0;
    B=3407.1*4.2/8.314;
    C=0;
    D=0;
    
elseif i==2 %parameter ji
        A=0;
        B=1662*4.2/8.314;
        C=0;
        D=0;
   
end

tij = A+B/T+C*log(T)+D*T;
   end

   %% critical temperature and composition
   
   function F = critVal(x,alpha)

% parameter for the substances

alpha_12=alpha;%NRTL(6,i);

tau_12 = tau(x(2),1); 
tau_21 = tau(x(2),2);
x1 = x(1);

% first and second derivative equals 0
% the first one is from henri el al
% the second via symbolyc derication in Matlab
F(1,1) = -1.0./(x1-1.0)+1.0./x1-tau_12.*exp(alpha_12.*tau_12.*-2.0).*1.0./(-x1+x1.*exp(-alpha_12.*tau_12)+1.0).^3.*2.0-tau_21.*exp(alpha_12.*tau_21.*-2.0).*1.0./(x1-exp(-alpha_12.*tau_21).*(x1-1.0)).^3.*2.0;
F(2,1) = 1.0./(x1-1.0).^2-1.0./x1.^2+tau_12.*exp(alpha_12.*tau_12.*-2.0).*(exp(-alpha_12.*tau_12)-1.0).*1.0./(-x1+x1.*exp(-alpha_12.*tau_12)+1.0).^4.*6.0-tau_21.*exp(alpha_12.*tau_21.*-2.0).*1.0./(x1-exp(-alpha_12.*tau_21).*(x1-1.0)).^4.*(exp(-alpha_12.*tau_21)-1.0).*6.0;
end

   %% LLE equilibrium
   function F=dge(x,T,alpha)

alpha_12 = alpha;%NRTL(6,i);


tau_12 = tau(T,1);
tau_21 = tau(T,2);

x_1=x(1,:);
x_1p=x(2,:);

x_2 = 1-x_1;
x_2p = 1-x_1p;

F(1,:) =  x_1.* exp(x_2.^2.*(tau_21.*(exp(-2.*alpha_12.*tau_21))./(x_1+x_2.*exp(-alpha_12.*tau_21)).^2+...
        tau_12.*(exp(-alpha_12.*tau_12))./(x_2+x_1.*exp(-alpha_12.*tau_12)).^2))-...
        x_1p.*exp(x_2p.^2.*(tau_21.*(exp(-2.*alpha_12.*tau_21))./(x_1p+x_2p.*exp(-alpha_12.*tau_21)).^2+...
        tau_12.*(exp(-alpha_12.*tau_12))./(x_2p+x_1p.*exp(-alpha_12.*tau_12)).^2));
    
F(2,:) = x_2.*exp(x_1.^2.*(tau_12.*(exp(-2.*alpha_12.*tau_12))./(x_2+x_1.*exp(-alpha_12.*tau_12)).^2+...
        tau_21.*(exp(-alpha_12.*tau_21))./(x_1+x_2.*exp(-alpha_12.*tau_21)).^2))-...
        x_2p.*exp(x_1p.^2.*(tau_12.*(exp(-2.*alpha_12.*tau_12))./(x_2p+x_1p.*exp(-alpha_12.*tau_12)).^2+...
        tau_21.*(exp(-alpha_12.*tau_21))./(x_1p+x_2p.*exp(-alpha_12.*tau_21)).^2));



   end

%% Target function

function t = target(x,T,alpha)

    z = dge(x,T,alpha);
    z1 = z(1,:);
    z2 = z(2,:);
    t = z1.^2 + z2.^2;

end