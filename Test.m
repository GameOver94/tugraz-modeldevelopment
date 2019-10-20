clc
clear
format long


%A Modified Sufficient Descent Polak–Ribiére–Polyak Type Conjugate Gradient Method for Unconstrained Optimization Problems

[x,y] = meshgrid(-5:0.5:5,-5:0.5:5);
z = func2({x,y});

xmin = min_cg(@func1, [3;3], 10^-3, 0.1, 50);


% figure()
% surf(x,y,z)
% title('surface')
%
figure()
hold on
contour(x,y,z)
plot(xmin(1,:),xmin(2,:),'-r*')
title('contour')
hold off


%% 

function z = func2(x)

% if ~ isvector(x)
%     error("Input for function is not a Vektor")
% end


z = 2*x{1}.^2 + 3*x{1}.*x{2} + 7.*x{2}.^2 + 8.*x{1} + 9.*x{2} + 10;
end

function z = func1(x)

% if ~ isvector(x)
%     error("Input for function is not a Vektor")
% end


z = 2*x(1,:).^2 + 3*x(1,:).*x(2,:) + 7.*x(2,:).^2 + 8.*x(1,:) + 9.*x(2,:) + 10;
end
%%

function [xmin] = min_cg(fun, x0, gtol, eps, maxiter)

k = 1;
x(:,k) = x0;



% Calculate Gradient at x0 y0

% generate matrix of function values for gradient calculation
xvec = ones(1, length(x0))*x(1,k);
for i = 2:length(x0)
    xvec = [xvec; ones(1, length(x0))*x(i,k)];
end

xvec_eps1 = xvec - eye(length(x0))*eps;
xvec_eps2 = xvec + eye(length(x0))*eps;


% calculate function values
fvec_eps1 = fun(xvec_eps1);
fvec_eps2 = fun(xvec_eps2);

% Gradient Vektor
g(:,k) = (fvec_eps2 - fvec_eps1)/(2*eps);


% initial search direction
d = -g;


while norm(g(:,k)) > gtol
    
    if k > maxiter
        disp('Maximum nubers of iteration reached')
        disp(['ABORDED! k= ', num2str(k)])
        break
    end

    % Step 2 :calculate alpha(k) with Armijo-Line search
    
    c = 0.5;
    tau = 0.9;
    
    m = g(:,k)'*d(:,k);
    
    t = -c*m;
    
    alpha(k) = 1;
    

    while fun(x(:,k))-fun(x(:,k)+alpha(k)*d(:,k)) < alpha(k)*t
        alpha(k) = tau * alpha(k);
    end

    
    % Step 3: calculate new x
    
    x(:,k+1) = x(:,k) + alpha(k)*d(:,k);
    
    % Step 4: calculate new serch direction
    
    mu= 0.5;
    
    % Calculate Gradient at xk+1 yk+1
%     df_dx(k+1) = diff([fun(x(k+1)-eps,y(k+1)),fun(x(k+1)+eps,y(k+1))])/(2 * eps);
%     df_dy(k+1) = diff([fun(x(k+1),y(k+1)-eps),fun(x(k+1),y(k+1)+eps)])/(2 * eps);
    
    % generate matrix of function values for gradient calculation
    xvec = ones(1, length(x0))*x(1,k+1);
    for i = 2:length(x0)
        xvec = [xvec; ones(1, length(x0))*x(i,k+1)];
    end
    
    xvec_eps1 = xvec - eye(length(x0))*eps;
    xvec_eps2 = xvec + eye(length(x0))*eps;
    
    
    % calculate function values
    fvec_eps1 = fun(xvec_eps1);
    fvec_eps2 = fun(xvec_eps2);
    
    % Gradient Vektor
    g(:,k+1) = (fvec_eps2 - fvec_eps1)/(2*eps);
    
    

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