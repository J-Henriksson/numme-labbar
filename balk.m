
clear; clc; close all;

H = 0.5;
g  = @(t) 8*exp(-t/2).*cos(3*t) - H;
gp = @(t) 8*exp(-t/2).*(-0.5*cos(3*t) - 3*sin(3*t));

tol   = 1e-8;
maxit = 100;

tN = newton(g, gp, 4.5, tol, maxit);
tS = secant(g, 4.2, 4.7, tol, maxit);

TN = tN(end);
TS = tS(end);

fprintf("Newton: T = %.12f, iterations = %d\n", TN, numel(tN)-1);
fprintf("Secant: T = %.12f, iterations = %d\n", TS, numel(tS)-2);

figure;
semilogy(abs(diff(tN)),'o-'); hold on;
semilogy(abs(diff(tS)),'s-');
grid on;
legend('Newton','Secant');
xlabel('n'); ylabel('|t_{n+1}-t_n|');

%Funktioner
function t_hist = newton(g, gp, t0, tol, maxit)
t_hist = t0;
for k = 1:maxit
    t = t_hist(end);
    t_new = t - g(t)/gp(t);
    t_hist(end+1) = t_new;
    if abs(t_new - t) < tol
        break;
    end
end
end

function t_hist = secant(g, t0, t1, tol, maxit)
t_hist = [t0 t1];
for k = 1:maxit
    a = t_hist(end-1); b = t_hist(end);
    denom = g(b) - g(a);
    if denom == 0
        error("Secant failed: division by zero.");
    end
    t_new = b - g(b)*(b-a)/denom;
    t_hist(end+1) = t_new; 
    if abs(t_new - b) < tol
        break;
    end
end
end
