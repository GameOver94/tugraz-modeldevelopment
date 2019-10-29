function [xmin] = min_cg(fun, x0, gtol, eps, maxiter, rho, delta, mu)
% Modified Polak-Rivére-Polyak


% Setp 0: initial Setup
k = 1;
x(:,k) = x0;

% Calculate Gradient at x0 y0

% generate matrix of function values for gradient calculation
xvec = ones(1, length(x0))*x(1,k);
for i = 2:length(x0)
    xvec = [xvec; ones(1, length(x0))*x(i,k)];
end

% xvec_eps1 = xvec - eye(length(x0))*eps;
% xvec_eps2 = xvec + eye(length(x0))*eps;
% 
% 
% % calculate function values
% fvec_eps1 = fun(xvec_eps1);
% fvec_eps2 = fun(xvec_eps2);
% 
% % Gradient Vektor
% g(:,k) = (fvec_eps2 - fvec_eps1)/(2*eps);

xvec_eps1 = xvec - 2.*eye(length(x0))*eps;
xvec_eps2 = xvec - eye(length(x0))*eps;
xvec_eps3 = xvec + eye(length(x0))*eps;
xvec_eps4 = xvec + 2.*eye(length(x0))*eps;


% calculate function values
fvec_eps1 = fun(xvec_eps1);
fvec_eps2 = fun(xvec_eps2);
fvec_eps3 = fun(xvec_eps3);
fvec_eps4 = fun(xvec_eps4);

% Gradient Vektor
g(:,k) = (fvec_eps1 - 8.*fvec_eps2 + 8.*fvec_eps3 - fvec_eps4) ./ (12*eps);

% initial search direction
d = -g;


while norm(g(:,k)) > gtol
    
    if k > maxiter
        disp('Maximum nubers of iteration reached')
        disp(['ABORDED! k= ', num2str(k)])
        break
    end
    
    % Step 2 :calculate alpha(k) with Armijo-Line search
    
    %rho = 0.5;
    %delta = 0.5;
    i = 1;
    
    alpha(k) = rho^i;       % Armijo
    %alpha(k) = 5;         % Backtracking
    
    while fun(x(:,k)+alpha(k)*d(:,k)) > fun(x(:,k))-delta*alpha(k)^2*norm(d(:,k))^2
        i=i+1;
        alpha(k)=rho^i;          % Armijo
        %alpha(k)=rho*alpha(k);   % Backtracking
    end
    
    
    % Step 3: calculate new x
    
    x(:,k+1) = x(:,k) + alpha(k)*d(:,k);
    
    % Step 4: calculate new serch direction
    
    %mu= 0.1;
      
    % generate matrix of function values for gradient calculation
    xvec = ones(1, length(x0))*x(1,k+1);
    for i = 2:length(x0)
        xvec = [xvec; ones(1, length(x0))*x(i,k+1)];
    end
    
%     xvec_eps1 = xvec - eye(length(x0))*eps;
%     xvec_eps2 = xvec + eye(length(x0))*eps;
%     
%     
%     % calculate function values
%     fvec_eps1 = fun(xvec_eps1);
%     fvec_eps2 = fun(xvec_eps2);
%     
%     % Gradient Vektor
%     g(:,k+1) = (fvec_eps2 - fvec_eps1)/(2*eps);


xvec_eps1 = xvec - 2.*eye(length(x0))*eps;
xvec_eps2 = xvec - eye(length(x0))*eps;
xvec_eps3 = xvec + eye(length(x0))*eps;
xvec_eps4 = xvec + 2.*eye(length(x0))*eps;


% calculate function values
fvec_eps1 = fun(xvec_eps1);
fvec_eps2 = fun(xvec_eps2);
fvec_eps3 = fun(xvec_eps3);
fvec_eps4 = fun(xvec_eps4);

% Gradient Vektor
g(:,k+1) = (fvec_eps1 - 8.*fvec_eps2 + 8.*fvec_eps3 - fvec_eps4) ./ (12*eps);
    
    
    
    if norm(g(:,k+1)) < gtol
        disp(['Minimum found after ', num2str(k), ' iterations.'])
        disp('The minnimum is at ')
        disp(x(:,k+1))
        break
    end
    
    % calculate beta
    beta(k+1) = g(:,k+1)'*(g(:,k+1)-g(:,k))/max(mu*norm(d(:,k))*norm(g(:,k+1)-g(:,k)),norm(g(:,k))^2);
    
    % calculate new search direction
    d(:,k+1) = -g(:,k+1)+beta(k+1)*d(:,k)-beta(k+1)*(g(:,k+1)'*d(:,k))/(g(:,k+1)'*(g(:,k+1)-g(:,k)))*(g(:,k+1)-g(:,k));
    
    % Step 5: Update k
    
    k = k+1;
    
end


xmin=x;

end

