% Model development
% WS 2019
%
% LLE MtOH - Hexane

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
c=fsolve(@(x)critVal(x,alpha),[0.1,283]);

T_crit = c(2);
x_crit = c(1);
T_min = -15+273.15;     % minimal valid temperature (-15°C bis 42°C)

%% LLE diagram

T(1)=T_crit;
x1_init(1) = x_crit-0.05; 
x1p_init(1) = x_crit+0.05; 


% Aspen dataset for comparing
A = xlsread('dataset_TvsMolFracMeOH(&hexane).xlsx');

x1Aspen = A(:,1);
T1Aspen = A(:,2);

% delete datasets which out of the temperature range
for i =1:length(T1Aspen)
    
    if T1Aspen(i,1)<=T_min-273.15
        x1Aspen(i) = nan;
        T1Aspen(i) = nan;
    end
end

s = 1

rho = 0.15;
for i=1:100
           
        %h=fsolve(@(x)dge(x,T(i),alpha),[x1_init(i);x1p_init(i)]);
        %min_cg(fun, x0, gtol, eps, maxiter, rho, delta, mu)
        h=min_cg_new(@(x)target(x,T(i),alpha),[x1_init(i);x1p_init(i)], 10^-7, 10^-5, 10^6, rho, 10^-2, 10^-1); % initial values !!
        x1(i) = h(1);
        x1p(i) = h(2);
        
        x1_init(i+1) = x1(i)-0.025;
        x1p_init(i+1) = x1p(i)+0.025;
        
        if s > 50
            disp('Update stepsize')
            rho = 10^-2.25;
        end
        
        
        temp=T;
         T(i+1) = T(i)-0.5;
         
         if T(i+1) < T_min
             break
         end       
         
         s=s+1
end
    

   plot(x_crit,T_crit-273.15,'go',x1,temp-273.15,'b--o',x1p,temp-273.15,'b--o',x1Aspen ,T1Aspen,'r--*')
   xlim([0 1]);
   title ('methanol-hexane');
   ylabel('temperature (°C)');
   xlabel('mol fraction');
   legend('critical Temperature','x1-left of T-crit','x1-right of T-crit','AspenSolution')

   
   %% tau calculation
   function tij = tau(T,i)

if i==1 %parameter ij
    A=-1.1544;
    B=734.514;
    C=0;
    D=0;
    
elseif i==2 %parameter ji
        A=-3.6511;
        B=1507.15;
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