%% Script ternary_RLLE_given_critKavalue
% Authors:  Group 6 - WS 2019/20
%           CHOWDHURY Hasan Mahmud - 11730325
%           HARRER Patrick - 01430527 
%           HIERZEGGER Robin - 01535430
%           KAIMBACHER Michael - 01431416 
%           SCHWINGSHACKL Julian - 01231490
% Purpose:  Estimate liquid-liquid equilibrium of A+B->C
%           with given Ka-value (-> calculate feed ratio A:B)

    clear variables; 
    clc;
    close all;

%% Calculate ternery LLE
% initial values
    x_1a(1) = 0.8;
    x_1b(1) = 0.2;
    x_3a(1) = 0;
    x_3b(1) = 0;
    K_a(1) = 0;
    T25 = 298.15;
    rho = 0.3;

    for i=1:40
        if i==1         
             h=amPRP(@(x)targetFunc(x,T25,K_a(i),0,0),[x_1a(i);x_1b(i);x_3a(i);x_3b(i)], ...
                10^-7, 10^6, rho, 10^-2, 10^-1);
        else
             h=amPRP(@(x)targetFunc(x,T25,K_a(i),0,0),[x_1a(i-1);x_1b(i-1);x_3a(i-1);x_3b(i-1)], ...
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
        K_a(i+1) = K_a(i) + 0.03;    
    end  
 
% Plot ternary diagram
%           A
%          /  \
%         /    \
%        B ---- C

    figure('Position', [50,50,1000,800]);
    ternplot(x_2a,x_3a,x_1a, 'r.-', 'majors', 10); 
    hold on
    ternplot(x_2b,x_3b,x_1b, 'b.-', 'majors', 10);  
    title('Reactive LLE: A + B <-> C', 'Position',[0.1 0.9])
    ternlabel('B','C','A'); 
    legend('Phase I', 'Phase II')
    
%% Calculate critical feed ratio - phi_F 
    x_1a_max = x_1a(1);
    x_2a_min = 1-x_1a_max;
    x_2b_max = x_2b(1);
    x_1b_min = 1-x_2b_max;

    phi_F_max = x_1a_max/x_2a_min;
    phi_F_min = x_1b_min/x_2b_max;
    
    clear('K_a')
    clear('x_1a','x_1b','x_2a','x_2b','x_3a','x_3b')

% given Ka-value
    K_a = 0.5;       %0.03<Ka_value<=0.75
    if K_a<=0.03
        error("Given K_a-value is too small because in terms of accuracy")
    elseif K_a>0.75
        error("K_a-value should be smaller than 0.76 because too close to the critical point")
    else
    end
    
% left branch = 0
% right branch = 1
    branch = [0,1];

% initial values
    x_1a = [x_1a_max,x_1a_max];
    x_1b = [x_1b_min,x_1b_min];
    x_3a = [0,0];
    x_3b = [0,0];
    
    phi_F = [phi_F_max,phi_F_min];
    phi = [x_1a_max,x_1b_min];
    x_1 = phi_F./(1+phi_F);

    rho = 0.4;

    for i=1:2   
         h=amPRP(@(x)targetFunc(x,T25,K_a,0,branch(i)),[x_1a(i);x_1b(i);x_3a(i);x_3b(i);phi_F(i);phi(i);x_1(i)],...
            10^-7, 10^6, rho, 10^-2, 10^-1);         
    
        x_1a(i) = h(1);
        x_1b(i) = h(2);
        x_3a(i) = h(3);
        x_3b(i) = h(4);   
        phi_F(i) = h(5);
        phi(i) = h(6);
        x_1(i) = h(7);    
        x_2a(i) = 1-x_1a(i)-x_3a(i);
        x_2b(i) = 1-x_1b(i)-x_3b(i);    
    end

    phi_F_crit = phi_F;

    if phi_F(2)<=phi_F_min || phi_F(1)>=phi_F_max     
        error("No miscibility gap for mixture");
    else
        display("Mixture has a Miscibility gap!")
    end
          
% tie lines
    X = [x_1a(1:3:end);x_1b(1:3:end)];
    Y = [x_2a(1:3:end);x_2b(1:3:end)];
    Z = [x_3a(1:3:end);x_3b(1:3:end)];
    
% plot ternary diagram 
    ternplot(Y,Z,X, 'go-');
    ternplot(x_1./phi_F_crit,(1-x_1-x_1./phi_F_crit),x_1,'s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
    ternplot(1./(1+phi_F_crit),[0,0],phi_F_crit./(1+phi_F_crit),'s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor',[.6 .6 1])
    legend('Phase I', 'Phase II', 'Tie Lines','K_{a,crit} (given)','feed ratio \phi_F (calculated)')
    [xpos,ypos] = terncoords(x_1(1)/phi_F_crit(1),(1-x_1(1)-x_1(1)/phi_F_crit(1)),x_1(1));
    text(xpos - 0.12,ypos + 0.015,['K_{a,crit} = ',num2str(K_a)],'Color','red','FontSize',9);
    [xpos,ypos] = terncoords(1./(1+phi_F_crit),[0,0],phi_F_crit./(1+phi_F_crit));
    text(xpos(1),ypos(1) - 0.035,['\phi_{F} = ',num2str(round(phi_F_crit(1),2))],'Color','red','FontSize',9);
    text(xpos(2),ypos(2) - 0.035,['\phi_{F} = ',num2str(round(phi_F_crit(2),2))],'Color','blue','FontSize',9);

%% Calculate progression of Reaction
    clear('K_a','phi_F')
    clear('x_1a','x_1b','x_2a','x_2b','x_3a','x_3b')
%% left branch
% initial values
    x_1a(1) = x_1a_max;
    x_1b(1) = x_1b_min;
    x_3a(1) = 0;
    x_3b(1) = 0;
    phi_F = phi_F_crit(1);
    phi(1) = x_1a_max;
    x_1(1) = phi_F/(1+phi_F);
    K_a(1) = 0;

    rho = 0.4;

    for i=1:40
        if i==1
            h=amPRP(@(x)targetFunc(x,T25,K_a(i),phi_F,0),[x_1a(i);x_1b(i);x_3a(i);x_3b(i);phi(i);x_1(i)],...
                10^-7, 10^6, rho, 10^-2, 10^-1);
        else
            h=amPRP(@(x)targetFunc(x,T25,K_a(i),phi_F,0),[x_1a(i-1);x_1b(i-1);x_3a(i-1);x_3b(i-1);phi(i-1);x_1(i-1)],...
                10^-7, 10^6, rho, 10^-2, 10^-1);
        end
         
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
        K_a(i+1) = K_a(i) + 0.03;          
    end

% plot line for reaction progress - left branch
    ternplot(x_1/phi_F,(1-x_1-x_1/phi_F),x_1,'m.-')

%% right branch
    clear('K_a','phi','x_1');
    clear('x_1a','x_1b','x_2a','x_2b','x_3a','x_3b')
% initial values
    x_1a(1) = x_1a_max;
    x_1b(1) = x_1b_min;
    x_3a(1) = 0;
    x_3b(1) = 0;
    phi_F = phi_F_crit(2);
    phi(1) = x_1b_min;
    x_1(1) = phi_F/(1+phi_F);
    K_a(1) = 0;

    rho = 0.4;

    for i=1:40
        if i==1
             h=amPRP(@(x)targetFunc(x,T25,K_a(i),phi_F,0),[x_1a(i);x_1b(i);x_3a(i);x_3b(i);phi(i);x_1(i)],...
                10^-7, 10^6, rho, 10^-2, 10^-1);
        else
             h=amPRP(@(x)targetFunc(x,T25,K_a(i),phi_F,0),[x_1a(i-1);x_1b(i-1);x_3a(i-1);x_3b(i-1);phi(i-1);x_1(i-1)],...
                10^-7, 10^6, rho, 10^-2, 10^-1);
        end
          
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
        K_a(i+1) = K_a(i) + 0.03;         
    end

% plot line for reaction progress - right branch
    ternplot(x_1/phi_F,(1-x_1-x_1/phi_F),x_1,'m.-')

% tie lines
    X = [x_1a(1:3:end);x_1b(1:3:end)];
    Y = [x_2a(1:3:end);x_2b(1:3:end)];
    Z = [x_3a(1:3:end);x_3b(1:3:end)];

% plot tie lines 
    ternplot(Y,Z,X, 'go-');

% plot Ka-value labeling
    [xpos, ypos] = terncoords(x_2b(1:3:end), x_3b(1:3:end), x_1b(1:3:end));
    Ka_lable = K_a(1:3:end);

    for i = 1:length(xpos)   
        text(xpos(i) + 0.01,ypos(i) + 0.005,num2str(Ka_lable(i)),'Color','green','FontSize',8);   
    end

    legend('Phase I', 'Phase II', 'Tie Lines','K_{a,crit} (given)','feed ratio \phi_F (calculated)','reaction progress')

%% Target function    
    function t = targetFunc(x,T,K_a,phi_F,branch)
        t = vecnorm(dge3(x,T,K_a,phi_F,branch)).^2;     % sum of Squares
    end