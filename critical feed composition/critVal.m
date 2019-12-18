function F = critVal(x,alpha_12)

tau_12 = tau(x(2,:),12); 
tau_21 = tau(x(2,:),21);
x1 = x(1,:);

% first and second derivative equals 0
% the first one is from henri el al
% the second via symbolyc derication in Matlab

F(1,:) = -1.0./(x1-1.0)+1.0./x1-tau_12.*exp(alpha_12.*tau_12.*-2.0).*1.0./(-x1+x1.*exp(-alpha_12.*tau_12)+1.0).^3.*2.0-tau_21.*exp(alpha_12.*tau_21.*-2.0).*1.0./(x1-exp(-alpha_12.*tau_21).*(x1-1.0)).^3.*2.0;
F(2,:) = 1.0./(x1-1.0).^2-1.0./x1.^2+tau_12.*exp(alpha_12.*tau_12.*-2.0).*(exp(-alpha_12.*tau_12)-1.0).*1.0./(-x1+x1.*exp(-alpha_12.*tau_12)+1.0).^4.*6.0-tau_21.*exp(alpha_12.*tau_21.*-2.0).*1.0./(x1-exp(-alpha_12.*tau_21).*(x1-1.0)).^4.*(exp(-alpha_12.*tau_21)-1.0).*6.0;
end
