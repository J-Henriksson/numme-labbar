clear; clc;

beta = 0.2;
a = 0; b = 20;

nList = [50 100 200 400];

nRef = 6400;
Vref = volume_trap(beta, a, b, nRef);

V = zeros(size(nList));
h = zeros(size(nList));
err = zeros(size(nList));

for k = 1:numel(nList)
    n = nList(k);
    h(k) = (b - a)/n;
    V(k) = volume_trap(beta, a, b, n);
    err(k) = abs(Vref - V(k));
end

% Observed order p from log(err) ~ p*log(h) + const
p = polyfit(log(h), log(err), 1);
order_est = p(1);

fprintf('--- Vattentorn (beta = %.3f) ---\n', beta);
fprintf('Reference (n = %d): Vref = %.10f\n\n', nRef, Vref);
fprintf('%8s %12s %18s\n', 'n', 'h', '|Vref - Vn|');
for k = 1:numel(nList)
    fprintf('%8d %12.6g %18.6e\n', nList(k), h(k), err(k));
end
fprintf('\nObserved order (slope in log-log): p ≈ %.3f\n', order_est);

% Also print the "best" volume among the tested nList
[~, idxBest] = max(nList);
fprintf('Best from n=%d: V = %.10f m^3\n', nList(idxBest), V(idxBest));

% -------- local function --------
function V = volume_trap(beta, a, b, n)
    % Composite trapezoidal rule
    x = linspace(a, b, n+1);
    y = (exp(beta*x) + 8) ./ (1 + (x/5).^3);
    f = pi * (y.^2);
    h = (b - a)/n;
    V = h * (0.5*f(1) + sum(f(2:end-1)) + 0.5*f(end));
end