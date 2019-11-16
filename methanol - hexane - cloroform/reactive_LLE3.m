% Model development
% WS 2019
%
%
%
% Reaktive LLE

% comp1: Hexane
% comp2: Methanol
% comp3: Chloroform


clear variables; 
clc;
close all;
format long;
options = optimoptions('fsolve','StepTolerance',1e-10);
%%
%% calculation LLE

% critical Temperatur and composition for the binary system 
alpha12 = 0.2;

%c=fsolve(@(x)critVal(x,alpha12),[0.1;283],options);
%min_cg(fun, x0, gtol, eps, maxiter, rho, delta, mu)
c=min_cg_new(@(x)target_crit(x,alpha12),[0.5;283],10^-8, 10^-6, 10^6, 0.2, 10^-2, 10^-1);

T_crit = c(2);
x_crit = c(1);




%% LLE diagram ternary system
x_3p(1) = 0; % one x has to be define in order to solve the equations

% initial values
x_1p(1) = 0.2;
x_1dp(1) = 0.75;
x_3dp(1) = 0;

T_calc = 293.15;
rho = 0.3;

for i=1:30
        if i==1   
        %h=fsolve(@(x)dge3(x,T_calc,x_3p(i)),[x_1p(i);x_1dp(i);x_3dp(i)],options); % starting values for x1', x1'',x3''
        h=min_cg_new(@(x)target(x,T_calc,x_3p(i)),[x_1p(i);x_1dp(i);x_3dp(i)], 10^-7, 10^-5, 10^6, rho, 10^-2, 10^-1);
        else
            %h=fsolve(@(x)dge3(x,T_calc,x_3p(i)),[x_1p(i-1);x_1dp(i-1);x_3dp(i-1)],options); % starting values for x1'',x2'',x3'
            h=min_cg_new(@(x)target(x,T_calc,x_3p(i)),[x_1p(i-1);x_1dp(i-1);x_3dp(i-1)], 10^-7, 10^-5, 10^6, rho, 10^-2, 10^-1);
        end
        
        x_1p(i) = h(1);
        x_1dp(i) = h(2);
        x_3dp(i)  = h(3) ;
        
        x_2p(i) = 1-x_1p(i)-x_3p(i);
        x_2dp(i) = 1-x_1dp(i)-x_3dp(i);

        
        x_3p(i+1) = x_3p(i)+0.01;
        


        if and(x_1p(i)==0, x_1dp(i)==0)
            break
        end
       
       
         
end

   x_3p(end) = [];
   
   %tie Lines
   X = [x_1p(1:3:end);x_1dp(1:3:end)];
   Y = [x_2p(1:3:end);x_2dp(1:3:end)];
   Z = [x_3p(1:3:end);x_3dp(1:3:end)];
   
   
   figure('Position', [50,50,800,800]);
   
    ternplot(x_1p,x_2p,x_3p, 'r.-', 'majors', 10); 
    hold on
    ternplot(x_1dp,x_2dp,x_3dp, 'b.-', 'majors', 10);
    
    ternplot(X,Y,Z, 'go-');
    
    title('Ternary Diagram from Hexane-Methanol-Chloroform at 25°C')
    ternlabel('hexane','menthol','chloroform'); 
    legend('Phase I', 'Phase II', 'Tie Lines')
   
%% Target function 
    
function t = target(x,T,x_3p)

    z = dge3(x,T,x_3p);
    z1 = z(1,:);
    z2 = z(2,:);
    z3 = z(3,:);
    t = z1.^2 + z2.^2 + z3.^2;

end

function t = target_crit(x,alpha_12)

    z = critVal(x,alpha_12);
    z1 = z(1,:);
    z2 = z(2,:);
    t = z1.^2 + z2.^2;

end  
  
   
   
  