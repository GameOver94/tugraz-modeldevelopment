% Model development
% WS 2019
%
%
%
% Reaktive LLE


clear variables; 
clc;
close all;

options = optimoptions('fsolve','StepTolerance',1e-10,'functionTolerance',1e-10);

%% calculation LLE

% critical Temperatur and composition
alpha12 = 0.2;

c=fsolve(@(x)critVal(x,alpha12),[0.5,283],options);

T_crit = c(2);
x_crit = c(1);

%% LLE diagram

T(1)=T_crit;
x1(1) = x_crit;
x1p(1) = x_crit;
for i=2:50
           
        T(i) = T(i-1)-1;
        h=fsolve(@(x)dge3(x,T(i),0),[0.1,0.9],options);
        x1(i) = h(1);
        x1p(i) = h(2);
           
end
    
   plot(x1,T-273.15,'b--o',x1p,T-273.15,'b--o') 
   grid on
   xlim([0 1]);
   title ('hexane - methanol');
   ylabel('temperature (°C)');
   xlabel('mol fraction hexan');

   
  