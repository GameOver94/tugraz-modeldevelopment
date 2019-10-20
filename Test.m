clc
clear
format long


%A Modified Sufficient Descent Polak–Ribiére–Polyak Type Conjugate Gradient Method for Unconstrained Optimization Problems

[x,y] = meshgrid(-5:0.5:5,-5:0.5:5);
z = func1(x,y);

[xmin,ymin] = min_cg(@func1, 3, 3, 10^-5, 0.1, 50);


% figure()
% surf(x,y,z)
% title('surface')
%
figure()
hold on
contour(x,y,z)
plot(xmin,ymin,'-r*')
title('contour')
hold off



function z = func1(x,y)

%[a,b,c,d,e,f] = args{:};
z = 2*x.^2 + 3*x.*y + 7.*y.^2 + 8.*x + 9.*y + 10;
end


function [xmin,ymin] = min_cg(fun, x0, y0, gtol, eps, maxiter)

k = 1;
x(k) = x0;
y(k) = y0;


% Calculate Gradient at x0 y0
df_dx(k) = diff([fun(x(k)-eps,y(k)),fun(x(k)+eps,y(k))])/(2 * eps);
df_dy(k) = diff([fun(x(k),y(k)-eps),fun(x(k),y(k)+eps)])/(2 * eps);

% Gradient Vektor
g{k} = [df_dx(k); df_dy(k)];

% initial search direction
d{k} = -g{k};


while norm(g{k}) > gtol
    
    if k > maxiter
        disp('Maximum nubers of iteration reached')
        disp(['ABORDED! k= ', num2str(k)])
        break
    end

    % Step 2 :calculate alpha(k) with Armijo-Line search
    
    c = 0.5;
    tau = 0.9;
    
    %m = g{k}'*(d{k}./norm(d{k}));
    m = g{k}'*d{k};
    
    t = -c*m;
    
    alpha(k) = 1;
    
    while fun(x(k),y(k)) - fun(x(k)+alpha(k)*d{k}(1), y(k)+alpha(k)*d{k}(2)) < alpha(k)*t
        alpha(k) = tau * alpha(k);
    end
    
    % Step 3: calculate new x
    
    x(k+1) = x(k) + alpha(k)*d{k}(1);
    y(k+1) = y(k) + alpha(k)*d{k}(2);
    
    % Step 4: calculate new serch direction
    
    mu= 0.5;
    
    % Calculate Gradient at xk+1 yk+1
    df_dx(k+1) = diff([fun(x(k+1)-eps,y(k+1)),fun(x(k+1)+eps,y(k+1))])/(2 * eps);
    df_dy(k+1) = diff([fun(x(k+1),y(k+1)-eps),fun(x(k+1),y(k+1)+eps)])/(2 * eps);
    
    % Gradient Vektor
    g{k+1} = [df_dx(k+1); df_dy(k+1)];
    

    if norm(g{k+1}) < gtol
        disp(['Minimum found after ', num2str(k), ' iterations.'])
        disp(['The minnimum is at ', num2str(x(k+1)),' ', num2str(y(k+1))])
        break
    end
    
    % calculate beta
    beta(k+1) = g{k+1}'*(g{k+1}-g{k})/max(mu*norm(d{k})*norm(g{k+1}-g{k}),norm(g{k})^2);
    
    % calculate new search direction
    d{k+1} = -g{k+1}+beta(k+1)*d{k}-beta(k+1)*(g{k+1}'*d{k})/(g{k+1}'*(g{k+1}-g{k}))*(g{k+1}-g{k});
    
    % Step 5: Update k
    
    k = k+1;
    
end


xmin=x;
ymin=y;

end