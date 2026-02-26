clear; clc; close all;

f = @(t,y) sin(3*t) - 2*y;             
y_exact = @(t) (93/65)*exp(-2*t) - (3/13)*cos(3*t) + (2/13)*sin(3*t);

y0 = 1.2;

%% (a) kollar residual r(t)=y'(t)-f(t,y(t))
tcheck = linspace(0,8,2001);
yex = y_exact(tcheck);


yex_p = -(186/65)*exp(-2*tcheck) + (9/13)*sin(3*tcheck) + (6/13)*cos(3*tcheck);

res = yex_p - (sin(3*tcheck) - 2*yex);
fprintf('Max residual |y'' - (sin(3t)-2y)| on [0,8]: %.3e\n\n', max(abs(res)));


forwardEuler  = @(T,n) fe_forward(f, y0, T, n);
backwardEuler = @(T,n) fe_backward(y0, T, n);

%% (b) Forward Euler på [0,8]
T = 8;
nList = [50 100 200 400];

figure;
hold on;
for n = nList
    [t,y] = forwardEuler(T,n);
    plot(t,y,'-');
end
plot(tcheck, y_exact(tcheck), 'k', 'LineWidth', 1.5);
hold off;
grid on;
title('Forward Euler on [0,8] for n=50,100,200,400');
xlabel('t'); ylabel('y(t)');
legend('n=50','n=100','n=200','n=400','exact','Location','best');

%% (c) Forward Euler error at t=8 vs h, loglog
hList = T ./ nList;
errF = zeros(size(nList));
for k = 1:numel(nList)
    n = nList(k);
    [t,y] = forwardEuler(T,n);
    errF(k) = abs(y_exact(T) - y(end));
end

figure;
loglog(hList, errF, 'o-');
grid on;
title('Forward Euler: |y(8)-y_h(8)| vs h (loglog)');
xlabel('h'); ylabel('error');

pF = polyfit(log(hList), log(errF), 1);
fprintf('Forward Euler estimated order p ≈ %.3f (slope in loglog)\n', pF(1));

%% (d) Backward Euler on [0,8]
figure;
hold on;
for n = nList
    [t,y] = backwardEuler(T,n);
    plot(t,y,'-');
end
plot(tcheck, y_exact(tcheck), 'k', 'LineWidth', 1.5);
hold off;
grid on;
title('Backward Euler on [0,8] for n=50,100,200,400 + exact');
xlabel('t'); ylabel('y(t)');
legend('n=50','n=100','n=200','n=400','exact','Location','best');

% Backward Euler error and order
errB = zeros(size(nList));
for k = 1:numel(nList)
    n = nList(k);
    [t,y] = backwardEuler(T,n);
    errB(k) = abs(y_exact(T) - y(end));
end

figure;
loglog(hList, errB, 'o-');
grid on;
title('Backward Euler: |y(8)-y_h(8)| vs h (loglog)');
xlabel('h'); ylabel('error');

pB = polyfit(log(hList), log(errB), 1);
fprintf('Backward Euler estimated order p ≈ %.3f (slope in loglog)\n\n', pB(1));

%% (e) Stability test: Forward Euler on [0,80]
Tlong = 80;
nListLong = [50 100 400 800];

figure;
for k = 1:numel(nListLong)
    n = nListLong(k);
    [t,y] = forwardEuler(Tlong,n);
    subplot(2,2,k);
    plot(t,y,'-'); hold on;
    plot(t, y_exact(t), 'k', 'LineWidth', 1);
    hold off; grid on;
    title(sprintf('Forward Euler, T=80, n=%d (h=%.3g)', n, Tlong/n));
    xlabel('t'); ylabel('y');
end

%% (f) Stability test: Backward Euler on [0,80]
figure;
for k = 1:numel(nListLong)
    n = nListLong(k);
    [t,y] = backwardEuler(Tlong,n);
    subplot(2,2,k);
    plot(t,y,'-'); hold on;
    plot(t, y_exact(t), 'k', 'LineWidth', 1);
    hold off; grid on;
    title(sprintf('Backward Euler, T=80, n=%d (h=%.3g)', n, Tlong/n));
    xlabel('t'); ylabel('y');
end


function [t,y] = fe_forward(f, y0, T, n)
    h = T/n;
    t = linspace(0,T,n+1).';
    y = zeros(n+1,1);
    y(1) = y0;
    for i = 1:n
        y(i+1) = y(i) + h*f(t(i), y(i));
    end
end

function [t,y] = fe_backward(y0, T, n)
    % For this specific ODE: y' = sin(3t) - 2y
    % Backward Euler step can be solved explicitly:
    % y_{i+1} = (y_i + h*sin(3*t_{i+1})) / (1 + 2h)
    h = T/n;
    t = linspace(0,T,n+1).';
    y = zeros(n+1,1);
    y(1) = y0;
    for i = 1:n
        tp = t(i+1);
        y(i+1) = (y(i) + h*sin(3*tp)) / (1 + 2*h);
    end
end