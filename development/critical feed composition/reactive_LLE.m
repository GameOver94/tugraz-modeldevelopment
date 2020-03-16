% Model development
% WS 2019
%
% LLE MtOH - Hexane

clear variables; 
clc;
close all;

options = optimoptions('fsolve','StepTolerance',1e-10,'functionTolerance',1e-10);

%% calculation LLE

% critical Temperatur and composition
alpha12 = 0.2;

%c=fsolve(@(x)critVal(x,alpha12),[0.5,283],options);
%min_cg(fun, x0, gtol, eps, maxiter, rho, delta, mu)
c=min_cg_new(@(x)target_crit(x,alpha12),[0.5;283],10^-8, 10^-6, 10^6, 0.2, 10^-2, 10^-1);

T_crit = c(2);
x_crit = c(1);

%% LLE diagram

T_min = -15+273.15;     % minimal valid temperature (-15°C bis 42°C)

T(1)=T_crit;
x1(1) = x_crit;
x1p(1) = x_crit;

x1_init = x_crit -0.1;
x1p_init = x_crit +0.1;

rho = 0.15;
for i=2:50
        
        
        T(i) = T(i-1)-1;
        disp(['Step: ', num2str(i), ' Temp: ', num2str(T(i)+273.15)]);
        
              
        if T(i) < 273.5 + 5
            rho = 0.05;
        end
        %h=fsolve(@(x)dge3(x,T(i),0),[0.1;0.9],options);
        %min_cg(fun, x0, gtol, eps, maxiter, rho, delta, mu)
        h=min_cg_new(@(x)target(x,T(i),0),[x1_init;x1p_init], 10^-7, 10^-5, 10^6, rho, 10^-2, 10^-1);
        x1(i) = h(1);
        x1p(i) = h(2);
        
        x1_init = x1(i)-0.025;
        x1p_init = x1p(i)+0.025;
        
           
end
    
   plot(x1,T-273.15,'b--o',x1p,T-273.15,'b--o') 
   grid on
   xlim([0 1]);
   title ('hexane - methanol');
   ylabel('temperature (°C)');
   xlabel('mol fraction hexan');
   
   
%% Target function

function t = target(x,T,x_3p)

    z = dge3(x,T,x_3p);
    z1 = z(1,:);
    z2 = z(2,:);
    t = z1.^2 + z2.^2;

end

function t = target_crit(x,alpha_12)

    z = critVal(x,alpha_12);
    z1 = z(1,:);
    z2 = z(2,:);
    t = z1.^2 + z2.^2;

end

   
  