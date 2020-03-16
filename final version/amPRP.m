%% Function amPRP
% Authors:  Group 6 - WS 2019/20
%           CHOWDHURY Hasan Mahmud - 11730325
%           HARRER Patrick - 01430527 
%           HIERZEGGER Robin - 01535430
%           KAIMBACHER Michael - 01431416 
%           SCHWINGSHACKL Julian - 01231490
% Purpose:  Modified Polak-Ribére-Polyak conjugate gradient algorithm (mPRP)
%           to solve Equations numerical

 function [xmin] = amPRP(fun, x0, eps, kmax, rho, delta, mu)

% Modified Polak-Rivére-Polyak

% Setp 0: initial setup
    k = 1;
    x(:,k) = x0;
    
    % generate matrix of function values for gradient calculation
        xvec = ones(1, length(x0))*x(1,k);
        for i = 2:length(x0)
            xvec = [xvec; ones(1, length(x0))*x(i,k)];
        end

    % Symetric Gradient
        xvec_eps1 = xvec - eye(length(x0))*eps;
    xvec_eps2 = xvec + eye(length(x0))*eps;
 
    % calculate function values
        fvec_eps1 = fun(xvec_eps1);
        fvec_eps2 = fun(xvec_eps2); 
        g(:,k) = (fvec_eps2 - fvec_eps1)/(2*eps);

    % initial search direction
        d = -g;

%% Going into loop
    
% Step 1: Check whether the gradient is big enough for the while-loop
    while norm(g(:,k)) > eps
 
      % Check whether number of maximum iterations have been reached
        if k > kmax
            disp('Maximum numbers of iteration reached')
            disp(['ABORDED! k = ', num2str(k)])
            disp(['Accuracy = ', num2str(norm(g(:,k)))])
            disp(' ')
            xmin = x(:,k);
            break
        end   
        
% Step 2 :calculate alpha(k) with Armijo-Line search
    	i = 1;   
        alpha(k) = rho^i;   
        while fun(x(:,k)+alpha(k)*d(:,k)) > fun(x(:,k))-delta*alpha(k)^2*norm(d(:,k))^2
            i = i+1;
            alpha(k) = rho^i;  
        end   
        
% Step 3: calculate new x    
        x(:,k+1) = x(:,k) + alpha(k)*d(:,k); 
        
         if norm(x(:,k+1) - x(:,k)) < eps^2       
             disp("Step doesn't improve result significantly")      
         end
       
% Step 4: calculate new serch direction
         
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
        if norm(g(:,k+1)) <= eps
            disp(['Minimum found after ', num2str(k), ' iterations.'])
            disp('The minimum is at ')
            disp(x(:,k+1))
            xmin=x(:,k+1);
            break
        end    
    % calculate beta
        beta(k+1) = g(:,k+1)'*(g(:,k+1)-g(:,k))/max(mu*norm(d(:,k))*norm(g(:,k+1)-g(:,k)),norm(g(:,k))^2);  
        
    % calculate new search direction
        d(:,k+1) = -g(:,k+1)+beta(k+1)*d(:,k)-beta(k+1)*(g(:,k+1)'*d(:,k))/(g(:,k+1)'*(g(:,k+1)-g(:,k)))*(g(:,k+1)-g(:,k));    

% Step 5: Update k   
        k = k+1;
        
    end
   
        if isnan(xmin)
            disp('ERROR NAN')
        end
 end
