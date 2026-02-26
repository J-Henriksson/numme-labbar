clear; clc;

a = 0; b = 20;
Vtarget = 1500;

nInt = 4000;

% Initial guesses for secant 
beta0 = 0.10;
beta1 = 0.30;

tol = 1e-10;
maxIter = 50;

F0 = F(beta0, a, b, nInt, Vtarget);
F1 = F(beta1, a, b, nInt, Vtarget);

fprintf('Secant method for V(beta) = %.1f\n', Vtarget);
fprintf('Using trapezoid with n = %d\n\n', nInt);
fprintf('%4s %12s %18s %18s\n', 'k', 'beta', 'F(beta)', '|beta_k - beta_{k-1}|');

fprintf('%4d %12.8f %18.8e %18s\n', 0, beta0, F0, '-');
fprintf('%4d %12.8f %18.8e %18.8e\n', 1, beta1, F1, abs(beta1 - beta0));

for k = 2:maxIter
    denom = (F1 - F0);
    if abs(denom) < 1e-16
        error('Secant failed: F(beta_k) - F(beta_{k-1}) too small.');
    end

    beta2 = beta1 - F1*(beta1 - beta0)/denom;
    F2 = F(beta2, a, b, nInt, Vtarget);

    fprintf('%4d %12.8f %18.8e %18.8e\n', k, beta2, F2, abs(beta2 - beta1));

    if abs(beta2 - beta1) < tol
        beta_star = beta2;
        V_star = volume_trap(beta_star, a, b, nInt);
        fprintf('\nConverged!\n');
        fprintf('beta ≈ %.10f\n', beta_star);
        fprintf('V(beta) ≈ %.10f m^3\n', V_star);
        return;
    end

    beta0 = beta1; F0 = F1;
    beta1 = beta2; F1 = F2;
end

warning('Reached max iterations without meeting tolerance.');
beta_star = beta1;
V_star = volume_trap(beta_star, a, b, nInt);
fprintf('Last beta ≈ %.10f, V ≈ %.10f\n', beta_star, V_star);

%Helper functions
function val = F(beta, a, b, n, Vtarget)
    val = volume_trap(beta, a, b, n) - Vtarget;
end

function V = volume_trap(beta, a, b, n)
    x = linspace(a, b, n+1);
    y = (exp(beta*x) + 8) ./ (1 + (x/5).^3);
    f = pi * (y.^2);
    h = (b - a)/n;
    V = h * (0.5*f(1) + sum(f(2:end-1)) + 0.5*f(end));
end