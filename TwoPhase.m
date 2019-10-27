clc
clear
format long


% Values for Contour Plot
[x,y] = meshgrid(-5:0.1:5,-5:0.1:5);
z = systemDraw({x,y});

fsolve(@systemSolve,[3;3])

% Calculate minnimum
%min_cg(fun, x0, gtol, eps, maxiter, rho, delta, mu)
%xmin = min_cg(@func2, [3;3], 10^-8, 0.1, 200, 0.5, 0.6, 0.1);


% Surface Plot
figure()
hold on
surf(x,y,z{1})
surf(x,y,z{2})
title('surface')
hold off

% Contour Plot
% figure()
% hold on
% contourf(x,y,z,15)
% plot(xmin(1,:),xmin(2,:),'-r*')
% title('contour')
% hold off


%% System of Equations Test Fucntion

function z = system(x)
    z(1)=4*x(1,:).^2+4*x(2,:).^2;
    z(2)=(x(1,:)-5.)^2+(x(2,:)-5).^2;
end

function z = systemSolve(x)
    z(1)=4*x(1).^2+4*x(2).^2;
    z(2)=(x(1)-5.)^2+(x(2)-5).^2;
end

function z = systemDraw(x)
    z{1}=4*x{1}.^2+4*x{2}.^2;
    z{2}=(x{1}-5).^2+(x{2}-5).^2;
end


%% Quadratic Test Function

function z = func1(x)
z = 2*x{1}.^2 + 3*x{1}.*x{2} + 7.*x{2}.^2 + 8.*x{1} + 9.*x{2} + 10;
end

function z = func2(x)
z = 2*x(1,:).^2 + 3*x(1,:).*x(2,:) + 7.*x(2,:).^2 + 8.*x(1,:) + 9.*x(2,:) + 10;
end