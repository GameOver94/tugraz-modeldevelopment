%% Script ternary_LLE_MetOH_Hex_Chloroform
% Authors:  Group 6 - WS 2019/20
%           CHOWDHURY Hasan Mahmud - 11730325
%           HARRER Patrick - 01430527 
%           HIERZEGGER Robin - 01535430
%           KAIMBACHER Michael - 01431416 
%           SCHWINGSHACKL Julian - 01231490
% Purpose:  Estimate liquid-liquid equilibrium  
%           of ternary system - n-hexane(1)/methanol(2)/chloroform(3)

    clear variables; 
    clc;
    close all;

%% Calculate ternery LLE
% initial values
    x_1a(1) = 0.8;
    x_1b(1) = 0.2;
    x_3a(1) = 0;
    x_3b(1) = 0;
    T25 = 298.15;   % 25°C
    rho = 0.3;
    K(1) = 0;       % as control variable

    for i=1:40
        if i==1         
            h=amPRP(@(x)targetFunc(x,T25,K(i)),[x_1a(i);x_1b(i);x_3a(i);x_3b(i)], ...
                10^-7, 10^6, rho, 10^-2, 10^-1);
        else
             h=amPRP(@(x)targetFunc(x,T25,K(i)),[x_1a(i-1);x_1b(i-1);x_3a(i-1);x_3b(i-1)], ...
                10^-7, 10^6, rho, 10^-2, 10^-1);
        end
       
        if abs(h(1)-h(2))<0.01
            disp(i)
            break
        end
    
        x_1a(i) = h(1);
        x_1b(i) = h(2);
        x_3a(i) = h(3);
        x_3b(i) = h(4);  
        x_2a(i) = 1-x_1a(i)-x_3a(i);
        x_2b(i) = 1-x_1b(i)-x_3b(i);   
        K(i+1) = K(i) + 0.03;    
    end  
   
   %tie Lines
   X = [x_1a(1:3:end);x_1b(1:3:end)];
   Y = [x_2a(1:3:end);x_2b(1:3:end)];
   Z = [x_3a(1:3:end);x_3b(1:3:end)];
 
   figure('Position', [50,50,1000,800]);
 
% Plot ternary diagram   
%       chloroform
%          /  \
%         /    \
% n-hexane ---- methanol
   
    ternplot(x_2a,x_3a,x_1a, 'r.-', 'majors', 10); 
    hold on
    ternplot(x_2b,x_3b,x_1b, 'b.-', 'majors', 10);    
    ternplot(Y,Z,X, 'go-');   
    title('Ternary Diagram from Hexane-Methanol-Chloroform at 25°C', 'Position',[0.1 0.9])
    ternlabel('methanol','chloroform','hexane'); 
    legend('Phase I', 'Phase II', 'Tie Lines')

%% Target function    
    function t = targetFunc(x,T,K) 
        t = vecnorm(dge3(x,T,K)).^2; % sum of Squares
    end