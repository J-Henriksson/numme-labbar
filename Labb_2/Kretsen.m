clear; clc; close all;

%start values
u0 = [1, 0];
tspan = [0, 20];

%constants for dampened circuit
L = 2;
C = 0.5;
R = 1;
[t, u] = ode45(@(t, u) odefnc(t, u, R, L, C), tspan, u0);

hold on
plot(t, u(:, 1), 'r')

%undampened circuit, (resitance = 0)
R = 0;
[t, u] = ode45(@(t, u) odefnc(t, u, R, L, C), tspan, u0);

figure(1)
plot(t, u(:,1), 'b')

title('lösning till ode:er')
xlabel('tid');
ylabel('laddning');
legend("R=1", "R=0")
hold off

% d) --------------------------------------------

tspan = [0, 40];
n = [40, 80, 160, 320];
R = 1;

figure(2)
[t1, u1] = ode45(@(t, u) odefnc(t, u, R, L, C), tspan, u0);
for i = 1:length(n)
    subplot(2,3,i);
    [t2, u2] = eulerfram(@(t, u) odefnc(t, u, R, L, C), tspan, n(i), u0);
    
    hold on
    plot(t1, u1(:, 1), 'r');
    plot(t2, u2(:, 1),'b');
    title(['n = ' num2str(n(i))])
    hold off
end


% ~ = t (not used in this function)
% u = derivative vector
% R, L, C = constants
function F = odefnc(~, u, R, L, C)
    F = zeros(2,1);
    F(1) = u(2);
    F(2) = -(R*u(2) + 1/C*u(1))/L;
end

% func = derivate function / vector of functions
% span = start/end values
% n = number of points
% u0 = start value(s)
function [t, u] = eulerfram(func, span, n, u0)
    h = (span(2) - span(1))/n;
    t = linspace(span(1), span(2), n+1);
    u = zeros(length(u0), n+1);
    u(:,1) = u0;
    for i = 1:n
    u(:,i+1) = u(:,i) + h*func(t(i), u(:,i));
    end
    u
    u = u'; %transpose to match ode45 ouput
end