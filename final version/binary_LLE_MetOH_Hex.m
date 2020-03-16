%% Script binary_LLE_MetOH_Hex
% Authors:  Group 6 - WS 2019/20
%           CHOWDHURY Hasan Mahmud - 11730325
%           HARRER Patrick - 01430527 
%           HIERZEGGER Robin - 01535430
%           KAIMBACHER Michael - 01431416 
%           SCHWINGSHACKL Julian - 01231490
% Purpose:  Estimate liquid-liquid equilibrium  
%           of binary system - methanol (1)/n-hexane (2)

    clear variables; 
    clc;
    close all;
    
%% Define parameters
    alpha = 0.2; %methanol-hexane (=Cij) alpha12 = alpha21
    
%% Calculate critical Temperatur and composition
    c = amPRP(@(x)targetFunc(x,0,alpha),[0.5;283],...
        10^-8, 10^6, 0.2, 10^-2, 10^-1);
    T_crit = c(2);
    x_crit = c(1);
    T_min = -15+273.15;     % minimal valid temperature(-15°C bis 42°C)
    
%% Calculate LLE
% phase 1:a
% phase 2:b

    T(1) = T_crit;
    x1a_init(1) = x_crit-0.05; 
    x1b_init(1) = x_crit+0.05; 
    rho = 0.3;
    s = 1;

    for i=1:100
        h=amPRP(@(x)targetFunc(x,T(i),alpha),[x1a_init(i);x1b_init(i)],...
            10^-5, 10^6, rho, 10^-2, 10^-1); 
        
        x1a(i) = h(1);
        x1b(i) = h(2);        
        x1a_init(i+1) = x1a(i)-0.025;
        x1b_init(i+1) = x1b(i)+0.025;
        
        if s > 50
            disp('Update stepsize')
            rho = 10^-2.25;
        end
        
        temp=T;
        T(i+1) = T(i)-0.5;
         
         if T(i+1) < T_min
             break
         end 
         s = s+1;
    end    

% plot diagram
   plot(x_crit,T_crit-273.15,'go',x1a,temp-273.15,'b--.',x1b,temp-273.15,'r--.')
   xlim([0 1]);
   title ('methanol-hexane');
   ylabel('temperature [°C]');
   xlabel('Methanol [mol-fraction]');
   legend('critical Temperature','x1-left of T-crit','x1-right of T-crit')
   
%% Target functions
    function t = targetFunc(x,T,alpha)
        t = vecnorm(dge2(x,T,alpha)).^2;    % sum of Squares
    end