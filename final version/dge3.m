%% Function dge3
% Authors:  Group 6 - WS 2019/20
%           CHOWDHURY Hasan Mahmud - 11730325
%           HARRER Patrick - 01430527 
%           HIERZEGGER Robin - 01535430
%           KAIMBACHER Michael - 01431416 
%           SCHWINGSHACKL Julian - 01231490
% Purpose:  EEstimate liquid-liquid equilibrium 
%           of ternary system - n-hexane(1)/methanol(2)/chloroform(3)
%           or Reactive LLE - of A+B->C

 function F=dge3(x,T,K,phi_F,branch)

% 1 - n-hexane
% 2 - methanol
% 3 - chloroform

    alpha_12 = 0.2;
    alpha_13 = 0.3;
    alpha_23 = 0.3;

    tau_12 = tau(T,12);
    tau_21 = tau(T,21);
    tau_13 = tau(T,13);
    tau_23 = tau(T,23);
    tau_31 = tau(T,31);
    tau_32 = tau(T,32);

% definition of the 3 variabes
    x_1a = x(1,:);
    x_1b = x(2,:);

    switch length(x)
        case 3 % Ternary LLE
            x_3b = 0;
        case 4 % Ternary LLE with reaction
            x_3a = x(3,:);
            x_3b  = x(4,:);
        case 6 % Reaction Progress
            x_3a = x(3,:);
            x_3b  = x(4,:);
            phi = x(5,:);
            x_1 = x(6,:);
        case 7 % Critical Ka and crit. feed ratio
            if (K == 0) && (phi_F ~= 0)
                x_3a = x(3,:);
                x_3b  = x(4,:);     
                K = x(5,:);
                phi = x(6,:);
                x_1 = x(7,:);
            elseif (K ~= 0) && (phi_F == 0)
                x_3a = x(3,:);
                x_3b  = x(4,:);
                phi_F = x(5,:);
                phi = x(6,:);
                x_1 = x(7,:);
            else
                error("No matching case found");
            end        
        otherwise
            error("No matching case found");
            return;
    end

% elimination of the other variables with standardization operation
    x_2a = 1-x_1a-x_3a;
    x_2b = 1-x_1b-x_3b;

% Parameters for NRTL Model 
    G12 = exp(-alpha_12*tau_12);
    G21 = exp(-alpha_12*tau_21);
    G13 = exp(-alpha_13*tau_13);
    G23 = exp(-alpha_23*tau_23);
    G31 = exp(-alpha_13*tau_31);
    G32 = exp(-alpha_23*tau_32);
    G11 = 1;
    G22 = 1;
    G33 = 1;
    tau_11 = 0;
    tau_22 = 0;
    tau_33 = 0;

% NRTL for activity cooeficients
    ln_gamma_1a = (x_1a.*tau_11.*G11+x_2a.*tau_21.*G21+x_3a.*tau_31.*G31)./(x_1a.*G11+x_2a.*G21+x_3a.*G31)+ ...
                  (x_1a.*G11)./(x_1a.*G11+x_2a.*G21+x_3a.*G31).*(tau_11-(x_1a.*tau_11.*G11+x_2a.*tau_21.*G21+x_3a.*tau_31.*G31)./(x_1a.*G11+x_2a.*G21+x_3a.*G31))+ ...
                  (x_2a.*G12)./(x_1a.*G12+x_2a.*G22+x_3a.*G32).*(tau_12-(x_1a.*tau_12.*G12+x_2a.*tau_22.*G22+x_3a.*tau_32.*G32)./(x_1a.*G12+x_2a.*G22+x_3a.*G32))+ ...
                  (x_3a.*G13)./(x_1a.*G13+x_2a.*G23+x_3a.*G33).*(tau_13-(x_1a.*tau_13.*G13+x_2a.*tau_23.*G23+x_3a.*tau_33.*G33)./(x_1a.*G13+x_2a.*G23+x_3a.*G33));

    ln_gamma_2a = (x_1a.*tau_12.*G12+x_2a.*tau_22.*G22+x_3a.*tau_32.*G32)./(x_1a.*G12+x_2a.*G22+x_3a.*G32)+ ...
                  (x_1a.*G21)./(x_1a.*G11+x_2a.*G21+x_3a.*G31).*(tau_21-(x_1a.*tau_11.*G11+x_2a.*tau_21.*G21+x_3a.*tau_31.*G31)./(x_1a.*G11+x_2a.*G21+x_3a.*G31))+ ...
                  (x_2a.*G22)./(x_1a.*G12+x_2a.*G22+x_3a.*G32).*(tau_22-(x_1a.*tau_12.*G12+x_2a.*tau_22.*G22+x_3a.*tau_32.*G32)./(x_1a.*G12+x_2a.*G22+x_3a.*G32))+ ...
                  (x_3a.*G23)./(x_1a.*G13+x_2a.*G23+x_3a.*G33).*(tau_23-(x_1a.*tau_13.*G13+x_2a.*tau_23.*G23+x_3a.*tau_33.*G33)./(x_1a.*G13+x_2a.*G23+x_3a.*G33));

    ln_gamma_3a = (x_1a.*tau_13.*G13+x_2a.*tau_23.*G23+x_3a.*tau_33.*G33)./(x_1a.*G13+x_2a.*G23+x_3a.*G33)+ ...
                  (x_1a.*G31)./(x_1a.*G11+x_2a.*G21+x_3a.*G31).*(tau_31-(x_1a.*tau_11.*G11+x_2a.*tau_21.*G21+x_3a.*tau_31.*G31)./(x_1a.*G11+x_2a.*G21+x_3a.*G31))+ ...
                  (x_2a.*G32)./(x_1a.*G12+x_2a.*G22+x_3a.*G32).*(tau_32-(x_1a.*tau_12.*G12+x_2a.*tau_22.*G22+x_3a.*tau_32.*G32)./(x_1a.*G12+x_2a.*G22+x_3a.*G32))+ ...
                  (x_3a.*G33)./(x_1a.*G13+x_2a.*G23+x_3a.*G33).*(tau_33-(x_1a.*tau_13.*G13+x_2a.*tau_23.*G23+x_3a.*tau_33.*G33)./(x_1a.*G13+x_2a.*G23+x_3a.*G33));

    ln_gamma_1b = (x_1b.*tau_11.*G11+x_2b.*tau_21.*G21+x_3b.*tau_31.*G31)./(x_1b.*G11+x_2b.*G21+x_3b.*G31)+ ...
                  (x_1b.*G11)./(x_1b.*G11+x_2b.*G21+x_3b.*G31).*(tau_11-(x_1b.*tau_11.*G11+x_2b.*tau_21.*G21+x_3b.*tau_31.*G31)./(x_1b.*G11+x_2b.*G21+x_3b.*G31))+ ...
                  (x_2b.*G12)./(x_1b.*G12+x_2b.*G22+x_3b.*G32).*(tau_12-(x_1b.*tau_12.*G12+x_2b.*tau_22.*G22+x_3b.*tau_32.*G32)./(x_1b.*G12+x_2b.*G22+x_3b.*G32))+ ...
                  (x_3b.*G13)./(x_1b.*G13+x_2b.*G23+x_3b.*G33).*(tau_13-(x_1b.*tau_13.*G13+x_2b.*tau_23.*G23+x_3b.*tau_33.*G33)./(x_1b.*G13+x_2b.*G23+x_3b.*G33));

    ln_gamma_2b = (x_1b.*tau_12.*G12+x_2b.*tau_22.*G22+x_3b.*tau_32.*G32)./(x_1b.*G12+x_2b.*G22+x_3b.*G32)+ ...
                  (x_1b.*G21)./(x_1b.*G11+x_2b.*G21+x_3b.*G31).*(tau_21-(x_1b.*tau_11.*G11+x_2b.*tau_21.*G21+x_3b.*tau_31.*G31)./(x_1b.*G11+x_2b.*G21+x_3b.*G31))+ ...
                  (x_2b.*G22)./(x_1b.*G12+x_2b.*G22+x_3b.*G32).*(tau_22-(x_1b.*tau_12.*G12+x_2b.*tau_22.*G22+x_3b.*tau_32.*G32)./(x_1b.*G12+x_2b.*G22+x_3b.*G32))+ ...
                  (x_3b.*G23)./(x_1b.*G13+x_2b.*G23+x_3b.*G33).*(tau_23-(x_1b.*tau_13.*G13+x_2b.*tau_23.*G23+x_3b.*tau_33.*G33)./(x_1b.*G13+x_2b.*G23+x_3b.*G33));

    ln_gamma_3b = (x_1b.*tau_13.*G13+x_2b.*tau_23.*G23+x_3b.*tau_33.*G33)./(x_1b.*G13+x_2b.*G23+x_3b.*G33)+ ...
                  (x_1b.*G31)./(x_1b.*G11+x_2b.*G21+x_3b.*G31).*(tau_31-(x_1b.*tau_11.*G11+x_2b.*tau_21.*G21+x_3b.*tau_31.*G31)./(x_1b.*G11+x_2b.*G21+x_3b.*G31))+ ...
                  (x_2b.*G32)./(x_1b.*G12+x_2b.*G22+x_3b.*G32).*(tau_32-(x_1b.*tau_12.*G12+x_2b.*tau_22.*G22+x_3b.*tau_32.*G32)./(x_1b.*G12+x_2b.*G22+x_3b.*G32))+ ...
                  (x_3b.*G33)./(x_1b.*G13+x_2b.*G23+x_3b.*G33).*(tau_33-(x_1b.*tau_13.*G13+x_2b.*tau_23.*G23+x_3b.*tau_33.*G33)./(x_1b.*G13+x_2b.*G23+x_3b.*G33));

    F(1,:) = x_1a.*exp(ln_gamma_1a)-x_1b.*exp(ln_gamma_1b);
    F(2,:) = x_2a.*exp(ln_gamma_2a)-x_2b.*exp(ln_gamma_2b);
    
    switch length(x)
        case 3
            F(3,:) = x_3a.*exp(ln_gamma_3a)-x_3b.*exp(ln_gamma_3b);
        
        case 4
            F(3,:) = x_3a.*exp(ln_gamma_3a)-x_3b.*exp(ln_gamma_3b);
            F(4,:) = x_3a.*exp(ln_gamma_3a)./(x_1a.*exp(ln_gamma_1a).*x_2a.*exp(ln_gamma_2a))-K;
        
        case 6
            F(3,:) = x_3a.*exp(ln_gamma_3a)-x_3b.*exp(ln_gamma_3b);
            F(4,:) = phi.*x_1a+(1-phi).*x_1b-x_1;
            F(5,:) = (phi.*x_2a+(1-phi).*x_2b)-x_1./phi_F;
            F(6,:) = x_3a.*exp(ln_gamma_3a)./(x_1a.*exp(ln_gamma_1a).*x_2a.*exp(ln_gamma_2a))-K;
        
        case 7
            if (branch == 0) % left branch
                F(3,:) = x_3a.*exp(ln_gamma_3a)-x_3b.*exp(ln_gamma_3b);
                F(4,:) = phi.*x_1a+(1-phi).*x_1b-x_1;
                F(5,:) = (phi.*x_2a+(1-phi).*x_2b)-x_1./phi_F;
                F(6,:) = x_3a.*exp(ln_gamma_3a)./(x_1a.*exp(ln_gamma_1a).*x_2a.*exp(ln_gamma_2a))-K;
                F(7,:) = phi-1;            
            elseif (branch == 1) %right branch
                F(3,:) = x_3a.*exp(ln_gamma_3a)-x_3b.*exp(ln_gamma_3b);
                F(4,:) = phi.*x_1a+(1-phi).*x_1b-x_1;
                F(5,:) = (phi.*x_2a+(1-phi).*x_2b)-x_1./phi_F;
                F(6,:) = x_3a.*exp(ln_gamma_3a)./(x_1a.*exp(ln_gamma_1a).*x_2a.*exp(ln_gamma_2a))-K;
                F(7,:) = phi;              
            else
                error("Branch index not valid");
            end       
    end
 end
