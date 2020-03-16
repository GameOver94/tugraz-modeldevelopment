%% Function dge2
% Authors:  Group 6 - WS 2019/20
%           CHOWDHURY Hasan Mahmud - 11730325
%           HARRER Patrick - 01430527 
%           HIERZEGGER Robin - 01535430
%           KAIMBACHER Michael - 01431416 
%           SCHWINGSHACKL Julian - 01231490
% Purpose:  Estimate liquid-liquid equilibrium 
%           of binary system - methanol (1)/n-hexane (2)

 function F=dge2(x,T,alpha)

 % parameter for the substances 
    alpha_12 = alpha; %alpha12 = alpha21

    if T == 0 % calculate critical point
        
            tau_12 = tau(x(2,:),12);
            tau_21 = tau(x(2,:),21);
            G12 = exp(-alpha_12.*tau_12);
            G21 = exp(-alpha_12.*tau_21);
            x1 = x(1,:);
            x2 = 1.-x1;
               
         % first derivative
            F(1,:) = (1./x1)+(1./x2)-(2.0.*G21.^2.*tau_21./((x1+x2.*G21).^2.*(x1+x2.*G21)))-...
                (2.0.*G12.^2.*tau_12./((x2+x1.*G12).^2.*(x2+x1.*G12)));    
         
         % second derivative
            F(2,:) = log((x1./x2))-((G21.^2.*tau_21)./((1-G21).*(x1+x2.*G21).^2))-...
               ((G12.^2.*tau_12)./((G12-1).*(x2+x1.*G12).^2));
                    
    elseif T ~= 0 % calculate binary LLE
           
            tau_12 = tau(T,12);
            tau_21 = tau(T,21);   
            G12 = exp(-alpha_12*tau_12);
            G21 = exp(-alpha_12*tau_21);
  
            x_1a = x(1,:);
            x_1b = x(2,:);
            x_2a = 1-x_1a;
            x_2b = 1-x_1b;
        
    % NRTL for activity cooeficients
         
    ln_gamma_1a = x_2a.^2.*((tau_21.*(G21./(x_1a+x_2a.*G21)).^2)+(tau_12.*G12./((x_2a+x_1a.*G12).^2)));
    ln_gamma_1b = x_2b.^2.*((tau_21.*(G21./(x_1b+x_2b.*G21)).^2)+(tau_12.*G12./((x_2b+x_1b.*G12).^2)));
    ln_gamma_2a = x_1a.^2.*((tau_12.*(G12./(x_2a+x_1a.*G12)).^2)+(tau_21.*G21./((x_1a+x_2a.*G21).^2)));
    ln_gamma_2b = x_1b.^2.*((tau_12.*(G12./(x_2b+x_1b.*G12)).^2)+(tau_21.*G21./((x_1b+x_2b.*G21).^2)));
    
       F(1,:) = x_1a.*exp(ln_gamma_1a)-x_1b.*exp(ln_gamma_1b);
       F(2,:) = x_2a.*exp(ln_gamma_2a)-x_2b.*exp(ln_gamma_2b);   
       
    else
                error("No matching case found");
            return;
    end   
 end