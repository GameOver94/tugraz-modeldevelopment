clc
clear
format long


%A Modified Sufficient Descent Polak�Ribi�re�Polyak Type Conjugate Gradient Method for Unconstrained Optimization Problems

[x,y] = meshgrid(-2:0.5:7,-2:0.5:7);
z = func2({x,y});

xmin = min_cg_debug(@func1, [5;5], 10^-3, 0.1, 50,0.25,0.5,0.3);


 figure()
 surf(x,y,z)
% title('surface')
%
figure()
hold on
contour(x,y,z)
for i=1:length(xmin(2,:))-1
    plot(xmin(1,i:i+1),xmin(2,i:i+1),'-r*','LineWidth',1.5)
    saveas(gcf,i + ".png")
    pause(0.5)
end
title('contour')
hold off


%% 

function z = func2(x)

% if ~ isvector(x)
%     error("Input for function is not a Vektor")
% end


%z = 2*x{1}.^2 + 3*x{1}.*x{2} + 7.*x{2}.^2 + 8.*x{1} + 9.*x{2} + 10;
% Booth function
z = (x{1}+2.*x{2}-7).^2+(2.*x{1}+x{2}-5).^2;

end

function z = func1(x)

% if ~ isvector(x)
%     error("Input for function is not a Vektor")
% end


%z = 2*x(1,:).^2 + 3*x(1,:).*x(2,:) + 7.*x(2,:).^2 + 8.*x(1,:) + 9.*x(2,:) + 10;
% Booth function
z = (x(1,:)+2.*x(2,:)-7).^2+(2.*x(1,:)+x(2,:)-5).^2;
end
% %%
% 
% function [xmin] = min_cg(fun, x0, gtol, eps, maxiter)
% 
% k = 1;
% x(:,k) = x0;
% 
% 
% 
% % Calculate Gradient at x0 y0
% 
% % generate matrix of function values for gradient calculation
% xvec = ones(1, length(x0))*x(1,k);
% for i = 2:length(x0)
%     xvec = [xvec; ones(1, length(x0))*x(i,k)];
% end
% 
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
% 
% 
% % initial search direction
% d = -g;
% 
% 
% while norm(g(:,k)) > gtol
%     
%     if k > maxiter
%         disp('Maximum nubers of iteration reached')
%         disp(['ABORDED! k= ', num2str(k)])
%         break
%     end
% 
%     % Step 2 :calculate alpha(k) with Standard Wolfe-Line Search Technique
%     
%     %c = 0.5;
%     tau = 0.9;
%     
%     %m = g(:,k)'*d(:,k);
%     
%     %t = -c*m;
%     
%     alpha(k) = 1;
%     
% 
%    % while fun(x(:,k))-fun(x(:,k)+alpha(k)*d(:,k)) < alpha(k)*t
%    %     alpha(k) = tau * alpha(k);
%    % end
% 
%    delta = 0.5;
%    delta1 = 0.25;
%    sigma = 0.8;
%    zeta = (delta*sigma/delta1)*0.5;
%    
%    % Pre Calcultion
%    v(:,k) = x(:,k) + alpha(k)*d(:,k);
%    e(:,k)=(-v(:,k)+x(:,k))/(alpha(k));
%    
%    
%     while fun(x(:,k)+alpha(k)*d(:,k))-fun(x(:,k)) > delta*alpha(k)*(g(:,k)'*d(:,k)) && e(:,k)'*d(:,k)<sigma*g(:,k)'*d(:,k)
%         alpha(k) = tau * alpha(k);
%         v(:,k) = x(:,k) + alpha(k)*d(:,k);
%         e(:,k)=(-v(:,k)+x(:,k))/(alpha(k));
%     end
%     
%     % Step 3: calculate new x
%     
%     %x(:,k+1) = x(:,k) + alpha(k)*d(:,k);
%     
%     
%     % New Calculation of:
%     v(:,k) = x(:,k) + alpha(k)*d(:,k);   
%     e(:,k)=(-v(:,k)+x(:,k))/(alpha(k));
%     
%     
%     if -delta1*(g(:,k)'*d(:,k)) >= delta*alpha(k)*(norm(d(:,k)))^2
%         x(:,k+1)= v(:,k);
%     else
%         x(:,k+1) = x(:,k)-(((e(:,k)'*(v(:,k)-x(:,k))+zeta*((norm(v(:,k)-x(:,k)))^2))/norm(g(:,k)))*g(:,k));
%         %x(:,k+1)= v(:,k);
%     end
%     
%     % Step 4: calculate new serch direction
%     
%     %mu= 0.5;
%     
%     % Calculate Gradient at xk+1 yk+1
% %     df_dx(k+1) = diff([fun(x(k+1)-eps,y(k+1)),fun(x(k+1)+eps,y(k+1))])/(2 * eps);
% %     df_dy(k+1) = diff([fun(x(k+1),y(k+1)-eps),fun(x(k+1),y(k+1)+eps)])/(2 * eps);
%     
%     % generate matrix of function values for gradient calculation
%     xvec = ones(1, length(x0))*x(1,k+1);
%     for i = 2:length(x0)
%         xvec = [xvec; ones(1, length(x0))*x(i,k+1)];
%     end
%     
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
%     
% 
%     if norm(g(:,k+1)) < gtol
%         disp(['Minimum found after ', num2str(k), ' iterations.'])
%         disp('The minimum is at ')
%         disp(x(:,k+1))
%         break
%     end
% 
%     
%     % calculate beta - modified version
%     %beta(k+1) = g(:,k+1)'*(g(:,k+1)-g(:,k))/max(mu*norm(d(:,k))*norm(g(:,k+1)-g(:,k)),norm(g(:,k))^2);
%     
%     beta(k+1) = g(:,k+1)'*(g(:,k+1)-g(:,k))/(norm(g(:,k+1)))^2; % Standard
%     
%     
%     % calculate new search direction - modified version
%     %d(:,k+1) = d(:,k+1) = -g(:,k+1)+beta(k+1)*d(:,k)-beta(k+1)*(g(:,k+1)'*d(:,k))/(g(:,k+1)'*(g(:,k+1)-g(:,k)))*(g(:,k+1)-g(:,k));
%  
%     d(:,k+1) = -g(:,k+1)+beta(k+1)*d(:,k);  % Standard 
%     
%     % Step 5: Update k
%     
%     k = k+1;
%     
% end
% 
% 
% xmin=x;
% 
% end