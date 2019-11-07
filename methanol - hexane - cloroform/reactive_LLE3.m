% Model development
% WS 2019
%
%
%
% Reaktive Extration

% comp1: Hexane
% comp2: Methanol
% comp3: Chloroform


clear variables; 
clc;
close all;
options = optimoptions('fsolve','StepTolerance',1e-10);
%%
%% calculation LLE

% critical Temperatur and composition for the binary system 
alpha12 = 0.2;

c=fsolve(@(x)critVal(x,alpha12),[0.1;283],options);


T_crit = c(2);
x_crit = c(1);




%% LLE diagram ternary system
x_3p(1) = 0; % one x has to be define in order to solve the equations

% initial values
x_1p(1) = 0.2;
x_1dp(1) = 0.75;
x_3dp(1) = 0;

T_calc = 293.15;

for i=1:31
        if i==1   
        
        h=fsolve(@(x)dge3(x,T_calc,x_3p(i)),[x_1p(i),x_1dp(i),x_3dp(i)],options); % starting values for x1', x1'',x3''
        else
            h=fsolve(@(x)dge3(x,T_calc,x_3p(i)),[x_1p(i-1),x_1dp(i-1),x_3dp(i-1)],options); % starting values for x1'',x2'',x3'
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
   X = [x_1p,x_1dp];
   Y = [x_2p,x_2dp];
   Z = [x_3p,x_3dp];
   
    ternplot(x_1p,x_2p,x_3p, 'r.', 'majors', 5); 
    hold on
    ternplot(x_1dp,x_2dp,x_3dp, 'b.', 'majors', 5);
    ternlabel('hexane','menthol','chloroform'); %heptane(1),toluene(2),sulfolane(3)
   
   
  
  
   
   
  