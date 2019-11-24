
function F=dge3(x,T,x_1p)

% x1' or one molefraction has to be specified

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
x_1dp = x(1,:);
x_3p = x(2,:);
x_3dp  = x(3,:);


% elimination of the other variables with standardization operation
x_2p = 1-x_1p-x_3p;
x_2dp = 1-x_1dp-x_3dp;


%NRTL Model for heptane (1), toluenee(2),sulfolane(3)
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
% theta_1 = (x1 + G21*x2 + G31*x3);
% theta_2 = (x2 + G12*x1 + G32*x3);
% theta_3 = (x3 + G13*x1 + G23*x2);
% 
% phi_1 = (x2*tau_21*G21 + x3*tau_31*G31);
% phi_2 = (x1*tau_12*G12 + x3*tau_32*G32);
% phi_3 = (x1*tau_13*G13 + x2*tau_23*G23);

%ln_gamma_1p = (x_1p + G21*x_2p + G31*x_3p)/(x_2p*tau_21*G21 + x_3p*tau_31*G31) - x_1p*(x_1p + G21*x_2p + G31*x_3p)/(x_2p*tau_21*G21 + x_3p*tau_31*G31)^2+x_2p*G12/(x_1p*tau_12*G12 + x_3p*tau_32*G32)*(tau_12-(x_2p + G12*x_1p + G32*x_3p)/(x_1p*tau_12*G12 + x_3p*tau_32*G32))+x_3p*G13/(x_1p*tau_13*G13 + x_2p*tau_23*G23)*(tau_13-(x_3p + G13*x_1p + G23*x_2p)/(x_1p*tau_13*G13 + x_2p*tau_23*G23));
%ln_gamma_2p = (x_2p + G12*x_1p + G32*x_3p)/(x_1p*tau_12*G12 + x_3p*tau_32*G32) - x_2p*(x_2p + G12*x_1p + G32*x_3p)/(x_1p*tau_12*G12 + x_3p*tau_32*G32)^2+x_1p*G21/(x_2p*tau_21*G21 + x_3p*tau_31*G31)*(tau_21-(x_1p + G21*x_2p + G31*x_3p)/(x_2p*tau_21*G21 + x_3p*tau_31*G31))+x_3p*G23/(x_1p*tau_13*G13 + x_2p*tau_23*G23)*(tau_23-(x_3p + G13*x_1p + G23*x_2p)/(x_1p*tau_13*G13 + x_2p*tau_23*G23));
%ln_gamma_3p = (x_3p + G13*x_1p + G23*x_2p)/(x_1p*tau_13*G13 + x_2p*tau_23*G23) - x_3p*(x_3p + G13*x_1p + G23*x_2p)/(x_1p*tau_13*G13 + x_2p*tau_23*G23)^2+x_1p*G31/(x_2p*tau_21*G21 + x_3p*tau_31*G31)*(tau_31-(x_1p + G21*x_2p + G31*x_3p)/(x_2p*tau_21*G21 + x_3p*tau_31*G31))+x_2p*G32/(x_1p*tau_12*G12 + x_3p*tau_32*G32)*(tau_32-(x_2p + G12*x_1p + G32*x_3p)/(x_1p*tau_12*G12 + x_3p*tau_32*G32));

%ln_gamma_1dp = (x_1dp + G21*x_2dp + G31*x_3dp)/(x_2dp*tau_21*G21 + x_3dp*tau_31*G31) - x_1dp*(x_1dp + G21*x_2dp + G31*x_3dp)/(x_2dp*tau_21*G21 + x_3dp*tau_31*G31)^2+x_2dp*G12/(x_1dp*tau_12*G12 + x_3dp*tau_32*G32)*(tau_12-(x_2dp + G12*x_1dp + G32*x_3dp)/(x_1dp*tau_12*G12 + x_3dp*tau_32*G32))+x_3dp*G13/(x_1dp*tau_13*G13 + x_2dp*tau_23*G23)*(tau_13-(x_3dp + G13*x_1dp + G23*x_2dp)/(x_1dp*tau_13*G13 + x_2dp*tau_23*G23));
%ln_gamma_2dp = (x_2dp + G12*x_1dp + G32*x_3dp)/(x_1dp*tau_12*G12 + x_3dp*tau_32*G32) - x_2dp*(x_2dp + G12*x_1dp + G32*x_3dp)/(x_1dp*tau_12*G12 + x_3dp*tau_32*G32)^2+x_1dp*G21/(x_2dp*tau_21*G21 + x_3dp*tau_31*G31)*(tau_21-(x_1dp + G21*x_2dp + G31*x_3dp)/(x_2dp*tau_21*G21 + x_3dp*tau_31*G31))+x_3dp*G23/(x_1dp*tau_13*G13 + x_2dp*tau_23*G23)*(tau_23-(x_3dp + G13*x_1dp + G23*x_2dp)/(x_1dp*tau_13*G13 + x_2dp*tau_23*G23));
%ln_gamma_3dp = (x_3dp + G13*x_1dp + G23*x_2dp)/(x_1dp*tau_13*G13 + x_2dp*tau_23*G23) - x_3dp*(x_3dp + G13*x_1dp + G23*x_2dp)/(x_1dp*tau_13*G13 + x_2dp*tau_23*G23)^2+x_1dp*G31/(x_2dp*tau_21*G21 + x_3dp*tau_31*G31)*(tau_31-(x_1dp + G21*x_2dp + G31*x_3dp)/(x_2dp*tau_21*G21 + x_3dp*tau_31*G31))+x_2dp*G32/(x_1dp*tau_12*G12 + x_3dp*tau_32*G32)*(tau_32-(x_2dp + G12*x_1dp + G32*x_3dp)/(x_1dp*tau_12*G12 + x_3dp*tau_32*G32));




ln_gamma_1p = (x_1p.*tau_11.*G11+x_2p.*tau_21.*G21+x_3p.*tau_31.*G31)./(x_1p.*G11+x_2p.*G21+x_3p.*G31)+ ...
              (x_1p.*G11)./(x_1p.*G11+x_2p.*G21+x_3p.*G31).*(tau_11-(x_1p.*tau_11.*G11+x_2p.*tau_21.*G21+x_3p.*tau_31.*G31)./(x_1p.*G11+x_2p.*G21+x_3p.*G31))+ ...
              (x_2p.*G12)./(x_1p.*G12+x_2p.*G22+x_3p.*G32).*(tau_12-(x_1p.*tau_12.*G12+x_2p.*tau_22.*G22+x_3p.*tau_32.*G32)./(x_1p.*G12+x_2p.*G22+x_3p.*G32))+ ...
              (x_3p.*G13)./(x_1p.*G13+x_2p.*G23+x_3p.*G33).*(tau_13-(x_1p.*tau_13.*G13+x_2p.*tau_23.*G23+x_3p.*tau_33.*G33)./(x_1p.*G13+x_2p.*G23+x_3p.*G33));
          
ln_gamma_2p = (x_1p.*tau_12.*G12+x_2p.*tau_22.*G22+x_3p.*tau_32.*G32)./(x_1p.*G12+x_2p.*G22+x_3p.*G32)+ ...
              (x_1p.*G21)./(x_1p.*G11+x_2p.*G21+x_3p.*G31).*(tau_21-(x_1p.*tau_11.*G11+x_2p.*tau_21.*G21+x_3p.*tau_31.*G31)./(x_1p.*G11+x_2p.*G21+x_3p.*G31))+ ...
              (x_2p.*G22)./(x_1p.*G12+x_2p.*G22+x_3p.*G32).*(tau_22-(x_1p.*tau_12.*G12+x_2p.*tau_22.*G22+x_3p.*tau_32.*G32)./(x_1p.*G12+x_2p.*G22+x_3p.*G32))+ ...
              (x_3p.*G23)./(x_1p.*G13+x_2p.*G23+x_3p.*G33).*(tau_23-(x_1p.*tau_13.*G13+x_2p.*tau_23.*G23+x_3p.*tau_33.*G33)./(x_1p.*G13+x_2p.*G23+x_3p.*G33));
          
ln_gamma_3p = (x_1p.*tau_13.*G13+x_2p.*tau_23.*G23+x_3p.*tau_33.*G33)./(x_1p.*G13+x_2p.*G23+x_3p.*G33)+ ...
              (x_1p.*G31)./(x_1p.*G11+x_2p.*G21+x_3p.*G31).*(tau_31-(x_1p.*tau_11.*G11+x_2p.*tau_21.*G21+x_3p.*tau_31.*G31)./(x_1p.*G11+x_2p.*G21+x_3p.*G31))+ ...
              (x_2p.*G32)./(x_1p.*G12+x_2p.*G22+x_3p.*G32).*(tau_32-(x_1p.*tau_12.*G12+x_2p.*tau_22.*G22+x_3p.*tau_32.*G32)./(x_1p.*G12+x_2p.*G22+x_3p.*G32))+ ...
              (x_3p.*G33)./(x_1p.*G13+x_2p.*G23+x_3p.*G33).*(tau_33-(x_1p.*tau_13.*G13+x_2p.*tau_23.*G23+x_3p.*tau_33.*G33)./(x_1p.*G13+x_2p.*G23+x_3p.*G33));
          

          
ln_gamma_1dp = (x_1dp.*tau_11.*G11+x_2dp.*tau_21.*G21+x_3dp.*tau_31.*G31)./(x_1dp.*G11+x_2dp.*G21+x_3dp.*G31)+ ...
              (x_1dp.*G11)./(x_1dp.*G11+x_2dp.*G21+x_3dp.*G31).*(tau_11-(x_1dp.*tau_11.*G11+x_2dp.*tau_21.*G21+x_3dp.*tau_31.*G31)./(x_1dp.*G11+x_2dp.*G21+x_3dp.*G31))+ ...
              (x_2dp.*G12)./(x_1dp.*G12+x_2dp.*G22+x_3dp.*G32).*(tau_12-(x_1dp.*tau_12.*G12+x_2dp.*tau_22.*G22+x_3dp.*tau_32.*G32)./(x_1dp.*G12+x_2dp.*G22+x_3dp.*G32))+ ...
              (x_3dp.*G13)./(x_1dp.*G13+x_2dp.*G23+x_3dp.*G33).*(tau_13-(x_1dp.*tau_13.*G13+x_2dp.*tau_23.*G23+x_3dp.*tau_33.*G33)./(x_1dp.*G13+x_2dp.*G23+x_3dp.*G33));
          
ln_gamma_2dp = (x_1dp.*tau_12.*G12+x_2dp.*tau_22.*G22+x_3dp.*tau_32.*G32)./(x_1dp.*G12+x_2dp.*G22+x_3dp.*G32)+ ...
              (x_1dp.*G21)./(x_1dp.*G11+x_2dp.*G21+x_3dp.*G31).*(tau_21-(x_1dp.*tau_11.*G11+x_2dp.*tau_21.*G21+x_3dp.*tau_31.*G31)./(x_1dp.*G11+x_2dp.*G21+x_3dp.*G31))+ ...
              (x_2dp.*G22)./(x_1dp.*G12+x_2dp.*G22+x_3dp.*G32).*(tau_22-(x_1dp.*tau_12.*G12+x_2dp.*tau_22.*G22+x_3dp.*tau_32.*G32)./(x_1dp.*G12+x_2dp.*G22+x_3dp.*G32))+ ...
              (x_3dp.*G23)./(x_1dp.*G13+x_2dp.*G23+x_3dp.*G33).*(tau_23-(x_1dp.*tau_13.*G13+x_2dp.*tau_23.*G23+x_3dp.*tau_33.*G33)./(x_1dp.*G13+x_2dp.*G23+x_3dp.*G33));
          
ln_gamma_3dp = (x_1dp.*tau_13.*G13+x_2dp.*tau_23.*G23+x_3dp.*tau_33.*G33)./(x_1dp.*G13+x_2dp.*G23+x_3dp.*G33)+ ...
              (x_1dp.*G31)./(x_1dp.*G11+x_2dp.*G21+x_3dp.*G31).*(tau_31-(x_1dp.*tau_11.*G11+x_2dp.*tau_21.*G21+x_3dp.*tau_31.*G31)./(x_1dp.*G11+x_2dp.*G21+x_3dp.*G31))+ ...
              (x_2dp.*G32)./(x_1dp.*G12+x_2dp.*G22+x_3dp.*G32).*(tau_32-(x_1dp.*tau_12.*G12+x_2dp.*tau_22.*G22+x_3dp.*tau_32.*G32)./(x_1dp.*G12+x_2dp.*G22+x_3dp.*G32))+ ...
              (x_3dp.*G33)./(x_1dp.*G13+x_2dp.*G23+x_3dp.*G33).*(tau_33-(x_1dp.*tau_13.*G13+x_2dp.*tau_23.*G23+x_3dp.*tau_33.*G33)./(x_1dp.*G13+x_2dp.*G23+x_3dp.*G33));

          
F(1,:) = x_1p.*exp(ln_gamma_1p)-x_1dp.*exp(ln_gamma_1dp);    
F(2,:) = x_2p.*exp(ln_gamma_2p)-x_2dp.*exp(ln_gamma_2dp);

if length(x) == 3
    %F(3,:) = x_3p.*exp(ln_gamma_3p)-x_3dp.*exp(ln_gamma_3dp);
    K = 4;
    %F(3,:) = K.*(x_1p.*exp(ln_gamma_1p).*x_2p.*exp(ln_gamma_2p))-(x_3dp.*exp(ln_gamma_3dp));
    F(3,:) = x_2p.*exp(ln_gamma_2p)./(K.*x_1p.*exp(ln_gamma_1p))-(x_3dp.*exp(ln_gamma_3dp));
end

end





  
   





