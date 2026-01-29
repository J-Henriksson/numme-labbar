clear; clc; close all;

%Settings
tol   = 1e-10;      % lämplig tolerans
maxit = 50;         % max amount of Newton iterations
showPlots = true;   % show circle plots + Newton iterates

%Data from assignment
A  = [175  950;
      410 2400;
      675 1730];

B  = [160 1008;
      381 2500;
      656 1760];

LA = [60; 75; 42];
LB = [45; 88; 57];

%Start guesses (near the right intersection for each pair of circles since that's the one we want).
%If Newton converges to the wrong intersection, adjust these slightly.
p0 = [205 1002;     % P1 start guess
      458 2458;     % P2 start guess
      712 1750];    % P3 start guess

%Newton
P = zeros(3,2);
iters = zeros(3,1);

for i = 1:3
    Ai = A(i,:);  Bi = B(i,:);
    rA = LA(i);   rB = LB(i);

    %Define the equation F([x;y]) = 0
    %F1 = distance^2 to A - LA^2
    %F2 = distance^2 to B - LB^2
    F = @(p) [ (p(1)-Ai(1))^2 + (p(2)-Ai(2))^2 - rA^2;
               (p(1)-Bi(1))^2 + (p(2)-Bi(2))^2 - rB^2 ];

    %Jacobian J(p) = dF/d[x,y]
    J = @(p) [ 2*(p(1)-Ai(1))  2*(p(2)-Ai(2));
               2*(p(1)-Bi(1))  2*(p(2)-Bi(2)) ];

    %Run Newton (no damping / no step-halving)
    [p_sol, p_hist, step_hist, res_hist] = newton2D(F, J, p0(i,:).', tol, maxit);

    P(i,:) = p_sol.';
    iters(i) = size(p_hist,2) - 1;

    %Print results
    fprintf('P%d:\n', i);
    fprintf('  start guess:   (%.6f, %.6f)\n', p0(i,1), p0(i,2));
    fprintf('  solution:      (%.12f, %.12f)\n', p_sol(1), p_sol(2));
    fprintf('  iterations:    %d\n', iters(i));
    fprintf('  final step:    %.3e\n', step_hist(end));
    fprintf('  final ||F||:   %.3e\n\n', res_hist(end));

    %Plot circles and Newton iterates (optional)
    if showPlots
        figure; hold on; axis equal; grid on;
        title(sprintf('P%d: circles + Newton iterates', i));
        xlabel('x'); ylabel('y');

        plotCircle(Ai, rA);
        plotCircle(Bi, rB);

        plot(Ai(1), Ai(2), 'k.', 'MarkerSize', 18);
        plot(Bi(1), Bi(2), 'k.', 'MarkerSize', 18);

        plot(p_hist(1,:), p_hist(2,:), 'o-', 'LineWidth', 1.5);
        plot(p_sol(1), p_sol(2), 'ks', 'MarkerSize', 10, 'LineWidth', 2);

        legend('Circle A','Circle B','A','B','Newton iterates','Solution', ...
               'Location','best');
    end
end

%Final coordinates
fprintf('=== Final coordinates ===\n');
fprintf('P1 = (%.12f, %.12f)\n', P(1,1), P(1,2));
fprintf('P2 = (%.12f, %.12f)\n', P(2,1), P(2,2));
fprintf('P3 = (%.12f, %.12f)\n', P(3,1), P(3,2));

%Interpolation (del b)
%Known points:
%P0 = (0,0)
%P4 = (1020,0)
P0 = [0 0];
P4 = [1020 0];

%Collect x and y coordinates for P0-P4
X = [P0(1), P(1,1), P(2,1), P(3,1), P4(1)];
Y = [P0(2), P(1,2), P(2,2), P(3,2), P4(2)];

%Sort points by x 
[X, idx] = sort(X);
Y = Y(idx);

%Find the 4th degree polynomial p(x) that goes through the 5 points
%p(x) = c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0
coeff = polyfit(X, Y, 4);

%Print coefficients
fprintf('\n=== 4th degree polynomial coefficients ===\n');
fprintf('p(x) = c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0\n');
fprintf('c4 = %.16e\n', coeff(1));
fprintf('c3 = %.16e\n', coeff(2));
fprintf('c2 = %.16e\n', coeff(3));
fprintf('c1 = %.16e\n', coeff(4));
fprintf('c0 = %.16e\n', coeff(5));

%Plot the road for x in (0,1020)
xx = linspace(0, 1020, 2000);
yy = polyval(coeff, xx);

figure; hold on; grid on;
plot(xx, yy, 'LineWidth', 1.5);      % the road p(x)
plot(X, Y, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);  % we mark interpolation points with 'o'
xlabel('x');
ylabel('y');
title('Road: 4th degree interpolation through P0..P4');
legend('p(x)', 'Points (o)', 'Location', 'best');

%Newton function (2D system)
function [p, hist, stepHist, resHist] = newton2D(F, J, p0, tol, maxit)
    p = p0;
    hist = p0;

    stepHist = [];
    resHist  = [];

    for k = 1:maxit
        Fp = F(p);
        Jp = J(p);

        %check if jacobian near singular
        if rcond(Jp) < 1e-14
            error('Jacobian near singular on iteration %d.', k);
        end

        %Newton step: solve J*delta = F
        delta = Jp \ Fp;

        %Update
        pNew = p - delta;

        %Store iteration history
        hist(:,end+1) = pNew;
        
        %step size and the residual
        step = norm(pNew - p, inf);
        res  = norm(F(pNew), 2);

        stepHist(end+1) = step;
        resHist(end+1)  = res;  

        %Stopping condition
        if step < tol
            p = pNew;
            return;
        end

        p = pNew;
    end
end

%Circle plotting function
function plotCircle(C, r)
    theta = linspace(0, 2*pi, 400);
    x = C(1) + r*cos(theta);
    y = C(2) + r*sin(theta);
    plot(x, y, 'LineWidth', 1.2);
end