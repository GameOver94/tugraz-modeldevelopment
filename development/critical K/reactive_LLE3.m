% Model development
% WS 2019
%
%
%
% Reaktive LLE

% comp1: A
% comp2: B
% comp3: C

% phase 1 : a
% pahse 2 : b

% Reaction: A+B-->C


clear variables; 
clc;
close all;
%format long;
options = optimoptions('fsolve','StepTolerance',1e-10);


%% calculation LLE

% critical Temperatur and composition for the binary system 
alpha12 = 0.2;

%c=fsolve(@(x)critVal(x,alpha12),[0.1;283],options);
%min_cg(fun, x0, gtol, eps, maxiter, rho, delta, mu)
c=min_cg_new(@(x)target_crit(x,alpha12),[0.5;283],10^-8, 10^-6, 10^6, 0.2, 10^-2, 10^-1);

T_crit = c(2);
x_crit = c(1);

%% Calculate ternery LLE

% initial values
x_1a(1) = 0.8;
x_1b(1) = 0.2;
x_3a(1) = 0;
x_3b(1) = 0;

K(1) = 0;

T_calc = 298.15;
rho = 0.3;

for i=1:40
    if i==1         
        %h=fsolve(@(x)dge3(x,T_calc,K(i)),[x_1a(i);x_1b(i);x_3a(i);x_3b(i),options);
        h=min_cg_new(@(x)target(x,T_calc,K(i),0),[x_1a(i);x_1b(i);x_3a(i);x_3b(i)], 10^-7, 10^-5, 10^6, rho, 10^-2, 10^-1);
    else
        %h=fsolve(@(x)dge3(x,T_calc,K(i)),[x_1a(i-1);x_1b(i-1);x_3a(i-1);x_3b(i-1);phi(i-1);x_1(i-1)],options);
        h=min_cg_new(@(x)target(x,T_calc,K(i),0),[x_1a(i-1);x_1b(i-1);x_3a(i-1);x_3b(i-1)], 10^-7, 10^-5, 10^6, rho, 10^-2, 10^-1);
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
   
   figure('Position', [50,50,1000,800]);
   
%       chloroform
%          /  \
%         /    \
%   hexane --- methanol
   
    ternplot(x_2a,x_3a,x_1a, 'r.-', 'majors', 10); 
    hold on
    ternplot(x_2b,x_3b,x_1b, 'b.-', 'majors', 10);
    
    title('ABC', 'Position',[0.1 0.9])
    ternlabel('B','C','A'); 
    legend('Phase I', 'Phase II')
  
    

%% Calculate Critical K

clear('K')
clear('x_1a','x_1b','x_2a','x_2b','x_3a','x_3b')

% ratio of the Feed
phi_F = 2;

% initial values
x_1a(1) = 0.8;
x_1b(1) = 0.2;
x_3a(1) = 0;
x_3b(1) = 0;

K(1) = 0;
phi(1) = 0.7;
x_1(1) = phi_F/(1+phi_F);


T_calc = 298.15;
rho = 0.5;

for i=1:1
    if i==1         
        %h=fsolve(@(x)dge3(x,T_calc,phi_F),[x_1a(i);x_1b(i);x_3a(i);x_3b(i);K(i);phi(i);x_1(i)],options); % starting values for x1', x1'',x3''
        h=min_cg_new(@(x)target(x,T_calc,0,phi_F),[x_1a(i);x_1b(i);x_3a(i);x_3b(i);K(i);phi(i);x_1(i)], 10^-7, 10^-6, 10^6, rho, 10^-2, 10^-1);
    else
        %h=fsolve(@(x)dge3(x,T_calc,phi_F),[x_1a(i-1);x_1b(i-1);x_3a(i-1);x_3b(i-1);K(i-1);phi(i-1);x_1(i-1)],options); % starting values for x1', x1'',x3''
        h=min_cg_new(@(x)target(x,T_calc,0,phi_F),[x_1a(i-1);x_1b(i-1);x_3a(i-1);x_3b(i-1);K(i-1);phi(i-1);x_1(i-1)], 10^-7, 10^-6, 10^6, rho, 10^-2, 10^-1);
    end
    
    
%     if abs(h(1)-h(2))<0.01
%         disp(i)
%         break
%     end
    
%     if or(h(5)<0,h(5)>1)
%         disp(i)
%         break
%     end
    
    x_1a(i) = h(1);
    x_1b(i) = h(2);
    x_3a(i) = h(3);
    x_3b(i) = h(4);
    
    K(i) = h(5);
    phi(i) = h(6);
    x_1(i) = h(7);

    x_2a(i) = 1-x_1a(i)-x_3a(i);
    x_2b(i) = 1-x_1b(i)-x_3b(i);
    
    
    
    
end

K_crit = K;
   
%tie Lines
X = [x_1a(1:3:end);x_1b(1:3:end)];
Y = [x_2a(1:3:end);x_2b(1:3:end)];
Z = [x_3a(1:3:end);x_3b(1:3:end)];

ternplot(Y,Z,X, 'go-');
ternplot(x_1/phi_F,(1-x_1-x_1/phi_F),x_1,'s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
legend('Phase I', 'Phase II', 'Tie Lines','crit. K')
[xpos,ypos] = terncoords(x_1/phi_F,(1-x_1-x_1/phi_F),x_1);
text(xpos - 0.15,ypos,['K = ',num2str(K_crit)],'Color','red','FontSize',12);


%% Calculate progression of Reaction

clear('K')
clear('x_1a','x_1b','x_2a','x_2b','x_3a','x_3b')

% initial values
x_1a(1) = 0.8;
x_1b(1) = 0.2;
x_3a(1) = 0;
x_3b(1) = 0;

phi(1) = 0.7;
x_1(1) = phi_F/(1+phi_F);


K(1) = 0;

T_calc = 298.15;
rho = 0.5;

for i=1:40
    if i==1
        %h=fsolve(@(x)dge3(x,T_calc,K(i),phi_F),[x_1a(i);x_1b(i);x_3a(i);x_3b(i);phi(i);x_1(i)],options); % starting values for x1', x1'',x3''
        h=min_cg_new(@(x)target(x,T_calc,K(i),phi_F),[x_1a(i);x_1b(i);x_3a(i);x_3b(i);phi(i);x_1(i)], 10^-7, 10^-5, 10^6, rho, 10^-2, 10^-1);
    else
        %h=fsolve(@(x)dge3(x,T_calc,K(i),phi_F),[x_1a(i-1);x_1b(i-1);x_3a(i-1);x_3b(i-1);phi(i-1);x_1(i-1)],options); % starting values for x1', x1'',x3''
        h=min_cg_new(@(x)target(x,T_calc,K(i),phi_F),[x_1a(i-1);x_1b(i-1);x_3a(i-1);x_3b(i-1);phi(i-1);x_1(i-1)], 10^-7, 10^-5, 10^6, rho, 10^-2, 10^-1);
    end
    
    
    %     if abs(h(1)-h(2))<0.01
    %         disp(i)
    %         break
    %     end
    
    if or(h(5)<0,h(5)>1)
        disp(i)
        break
    end
    
    x_1a(i) = h(1);
    x_1b(i) = h(2);
    x_3a(i) = h(3);
    x_3b(i) = h(4);
    
    phi(i) = h(5);
    x_1(i) = h(6);
    
    x_2a(i) = 1-x_1a(i)-x_3a(i);
    x_2b(i) = 1-x_1b(i)-x_3b(i);
    
    K(i+1) = K(i) + 0.03;
    
    
    
end


%tie Lines
X = [x_1a(1:3:end);x_1b(1:3:end)];
Y = [x_2a(1:3:end);x_2b(1:3:end)];
Z = [x_3a(1:3:end);x_3b(1:3:end)];

ternplot(x_1/phi_F,(1-x_1-x_1/phi_F),x_1,'b.-')
ternplot(Y,Z,X, 'go-');

legend('Phase I', 'Phase II', 'Tie Lines','crit. K','reaction progress')

%% Target function 
    
function t = target(x,T,K,phi_F)

z = dge3(x,T,K,phi_F);

if length(x) == 4
    z1 = z(1,:);
    z2 = z(2,:);
    z3 = z(3,:);
    z4 = z(4,:);
    t = z1.^2 + z2.^2 + z3.^2 + z4.^2;
elseif length(x) == 6
    z1 = z(1,:);
    z2 = z(2,:);
    z3 = z(3,:);
    z4 = z(4,:);
    z5 = z(5,:);
    z6 = z(6,:);
    t = z1.^2 + z2.^2 + z3.^2 + z4.^2 + z5.^2 + z6.^2;
elseif length(x) == 7
    z1 = z(1,:);
    z2 = z(2,:);
    z3 = z(3,:);
    z4 = z(4,:);
    z5 = z(5,:);
    z6 = z(6,:);
    z7 = z(7,:);
    t = z1.^2 + z2.^2 + z3.^2 + z4.^2 + z5.^2 + z6.^2 + z7.^2;
end

end

function t = target_crit(x,alpha_12)

    z = critVal(x,alpha_12);
    z1 = z(1,:);
    z2 = z(2,:);
    t = z1.^2 + z2.^2;

end  
  
   
   
  