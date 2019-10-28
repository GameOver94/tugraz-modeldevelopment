clc
clear
format long


% Values for Contour Plot
[x,y] = meshgrid(-5:0.5:5,-5:0.5:20);
z = systemDraw({x,y});

fsolve(@system,[3;3])
%min_cg(fun, x0, gtol, eps, maxiter, rho, delta, mu)
xmin = min_cg(@target, [4;10], 10^-8, 0.1, 200, 0.1, 0.5, 0.1);


% Calculate minnimum
%min_cg(fun, x0, gtol, eps, maxiter, rho, delta, mu)
%xmin = min_cg(@func2, [3;3], 10^-8, 0.1, 200, 0.5, 0.6, 0.1);


% Surface Plot
figure()
surf(x,y,z)
title('surface')


% Contour Plot
figure()
hold on
contourf(x,y,z,15)
plot(xmin(1,:),xmin(2,:),'-r*')
title('contour')
hold off


%% HIE82   https://www-m2.ma.tum.de/foswiki/pub/M2/Allgemeines/SemWs09/testprob.pdf

function z = system(x)
    x1 = x(1,:);
    x2 = x(2,:);
    z(1,:)=x2-10;
    z(2,:)=x1.*x2-5*10^-4;
end


function z = systemDraw(x)
    z1=x{2}-10;
    z2=x{1}.*x{2}-5*10^-4;
    
    z=sqrt(z1.^2+z2.^2);
end

function t = target(x)
    z = system(x);
    z1 = z(1,:);
    z2 = z(2,:);
    t = sqrt(z1.^2 + z2.^2);
end


%% Quadratic Test Function

function z = func1(x)
z = 2*x{1}.^2 + 3*x{1}.*x{2} + 7.*x{2}.^2 + 8.*x{1} + 9.*x{2} + 10;
end

function z = func2(x)
z = 2*x(1,:).^2 + 3*x(1,:).*x(2,:) + 7.*x(2,:).^2 + 8.*x(1,:) + 9.*x(2,:) + 10;
end