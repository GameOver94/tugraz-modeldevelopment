%% Script ternary_RLLE_calc_critKavalue
% Authors:  Group 6 - WS 2019/20
%           CHOWDHURY Hasan Mahmud - 11730325
%           HARRER Patrick - 01430527 
%           HIERZEGGER Robin - 01535430
%           KAIMBACHER Michael - 01431416 
%           SCHWINGSHACKL Julian - 01231490
% Purpose:  Estimate liquid-liquid equilibrium  of A+B->C
%           to calculate critical Ka-value (with given feed ratio A:B)

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
             h=amPRP(@(x)targetFunc(x,T25,K_a(i),0,0),[x_1a(i);x_1b(i);x_3a(i);x_3b(i)],...
                10^-7, 10^6, rho, 10^-2, 10^-1);
        else
             h=amPRP(@(x)targetFunc(x,T25,K_a(i),0,0),[x_1a(i-1);x_1b(i-1);x_3a(i-1);x_3b(i-1)],...
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
    
%% Calculate Critical Ka
    phi_F = 2; %ratio of the feed 0.2440 to 0.65 % 1.2 to 4.3377 
    x_1a_max = x_1a(1);
    x_2a_min = 1-x_1a_max;
    x_2b_max = x_2b(1);
    x_1b_min = 1-x_2b_max;

    phi_F_max = x_1a_max/x_2a_min;
    phi_F_min = x_1b_min/x_2b_max;
    
    if phi_F<=phi_F_min || phi_F>= phi_F_max   
        error("No miscibility gap for mixture");
    else
        display("Mixture has a Miscibility gap!")
    end
   
%Determine whether the distance to the critical point is large enough   
    s = length(x_1a);
%control which branch is affected
%set branch = 0 for left and branch = 1 for right branch
    if phi_F<(x_1b(s)/x_2b(s))-0.13
        branch = 1;      %right branch
        display("crit. Ka at right branch")
        phi_start = x_2a_min;
    elseif phi_F>(x_1a(s)/x_2a(s))+0.22
        branch = 0;      %left branch
        display("crit. Ka at left branch")
        phi_start = x_2b_max;
    else
        error("No solution, because too close to the critical point!");   
    end
 
% clear variables
    clear('K_a')
    clear('x_1a','x_1b','x_2a','x_2b','x_3a','x_3b')
        
% initial values
    x_1a(1) = x_1a_max;
    x_1b(1) = x_1b_min;
    x_3a(1) = 0;
    x_3b(1) = 0;    
    K_a(1) = 0;       
    x_1(1) = phi_F/(1+phi_F);
    phi(1) = phi_start;
    rho = 0.4;

    for i=1:1
        if i==1           
             h=amPRP(@(x)targetFunc(x,T25,0,phi_F,branch),[x_1a(i);x_1b(i);x_3a(i);x_3b(i);K_a(i);phi(i);x_1(i)],...
                10^-7, 10^6, rho, 10^-2, 10^-1); 
        else        
             h=amPRP(@(x)targetFunc(x,T25,0,phi_F,branch),[x_1a(i-1);x_1b(i-1);x_3a(i-1);x_3b(i-1);K_a(i-1);phi(i-1);x_1(i-1)],...
                10^-7, 10^6, rho, 10^-2, 10^-1); 
        end
        
        x_1a(i) = h(1);
        x_1b(i) = h(2);
        x_3a(i) = h(3);
        x_3b(i) = h(4);   
        K_a(i) = h(5);
        phi(i) = h(6);
        x_1(i) = h(7);
        x_2a(i) = 1-x_1a(i)-x_3a(i);
        x_2b(i) = 1-x_1b(i)-x_3b(i);   
    end
    
    Ka_crit = K_a;
   
% tie lines
    X = [x_1a(1:3:end);x_1b(1:3:end)];
    Y = [x_2a(1:3:end);x_2b(1:3:end)];
    Z = [x_3a(1:3:end);x_3b(1:3:end)];
    
% plot ternary diagram 
    ternplot(Y,Z,X, 'go-');
    ternplot(x_1/phi_F,(1-x_1-x_1/phi_F),x_1,'s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6]) 
    ternplot(1./(1+phi_F),[0,0],phi_F./(1+phi_F),'s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor',[.6 .6 1])
    legend('Phase I', 'Phase II', 'Tie Lines','K_{a,crit}','feed ratio \phi_F (given)')
    [xpos,ypos] = terncoords(x_1/phi_F,(1-x_1-x_1/phi_F),x_1);
    
    if branch == 0
    text(xpos - 0.12,ypos + 0.015,['K_{a,crit} = ',num2str(round(Ka_crit,2))],'Color','red','FontSize',9);   
    else 
    text(xpos + 0.01,ypos + 0.025,['K_{a,crit} = ',num2str(round(Ka_crit,2))],'Color','red','FontSize',9);   
    end
%% Calculate progression of Reaction
    clear('K_a')
    clear('x_1a','x_1b','x_2a','x_2b','x_3a','x_3b')
% initial values
    x_1a(1) = x_1a_max;
    x_1b(1) = x_1b_min;
    x_3a(1) = 0;
    x_3b(1) = 0;
    phi(1) = phi_start;
    x_1(1) = phi_F/(1+phi_F);
    K_a(1) = 0;
    rho = 0.4;

    for i=1:40
        if i==1      
             h=amPRP(@(x)targetFunc(x,T25,K_a(i),phi_F,branch),[x_1a(i);x_1b(i);x_3a(i);x_3b(i);phi(i);x_1(i)],...
                10^-7, 10^6, rho, 10^-2, 10^-1); 
        else
             h=amPRP(@(x)targetFunc(x,T25,K_a(i),phi_F,branch),[x_1a(i-1);x_1b(i-1);x_3a(i-1);x_3b(i-1);phi(i-1);x_1(i-1)],...
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

% tie lines
    X = [x_1a(1:3:end);x_1b(1:3:end)];
    Y = [x_2a(1:3:end);x_2b(1:3:end)];
    Z = [x_3a(1:3:end);x_3b(1:3:end)];

% plot in ternary diagram
    ternplot(x_1/phi_F,(1-x_1-x_1/phi_F),x_1,'b.-')
    ternplot(Y,Z,X, 'go-');
    legend('Phase I', 'Phase II', 'Tie Lines','K_{a,crit} (calculated)','feed ratio \phi_F (given)','reaction progress')   
    [xpos,ypos] = terncoords(x_1(1)/phi_F(1),-0.035,x_1(1));
    text(xpos(1),ypos(1)-0.005,['\phi_F = ',num2str(phi_F(1))],'Color','magenta','FontSize',9);
    
% plot K-value labeling
    [xpos, ypos] = terncoords(x_2b(1:3:end), x_3b(1:3:end), x_1b(1:3:end));
    Ka_lable = K_a(1:3:end);

    for i = 1:length(xpos)   
        text(xpos(i) + 0.01,ypos(i) + 0.005,num2str(Ka_lable(i)),'Color','green','FontSize',8);   
    end

    legend('Phase I', 'Phase II', 'Tie Lines','K_{a,crit} (calculated)','feed ratio \phi_F (given)','reaction progress')

%% Target function        
    function t = targetFunc(x,T,K_a,phi_F,branch)
        t = vecnorm(dge3(x,T,K_a,phi_F,branch)).^2;     % sum of Squares 
    end