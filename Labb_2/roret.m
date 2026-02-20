clear; clc; close all;

% Parametrar för Roret.m (Uppgift 4d med subplots)
Ti = 450;
Te = 20;
Ri = 1;
Ry = 2;
k = 1;

N = 10000;
h = (Ry - Ri) / (N + 1);

% r-vektorer för plottning och beräkning
r = Ri + h : h : Ry; 
r_full = [Ri, r]; % Inkluderar inre randen r=1 för plotten

% Vektor för alpha och för att spara resultaten
alpha_vec = linspace(0, 10, 21);
T_out = zeros(length(alpha_vec), 1); 

% --- BYGG DEN STATISKA DELEN AV SYSTEMET ---
A = sparse(N+1, N+1);
b = zeros(N+1, 1);

% Rad 1
A(1, 1) = -2 * r(1);
A(1, 2) = r(1) + h/2;
b(1)    = -(r(1) - h/2) * Ti;

% Rad 2 till N
for i = 2:N
    A(i, i-1) = r(i) - h/2;
    A(i, i)   = -2 * r(i);
    A(i, i+1) = r(i) + h/2;
    b(i)      = 0;
end

% --- FÖRBERED FIGUREN ---
% Skapar ett fönster som är lite bredare för att rymma två grafer
figure('Position', [100, 100, 1000, 400]); 

% Välj färgskala för profilerna (parula är MATLABs standard)
colors = parula(length(alpha_vec)); 

% Aktivera den högra subplotten för att rita profilerna i loopen
subplot(1, 2, 2);
hold on;

% --- LOOPA ÖVER OLIKA ALPHA ---
for m = 1:length(alpha_vec)
    alpha = alpha_vec(m);
    
    % Uppdatera sista raden
    A(N+1, N)   = 2 * r(N+1);
    A(N+1, N+1) = -2 * r(N+1) - (2*h*alpha/k) * (r(N+1) + h/2);
    b(N+1)      = -(2*h*alpha*Te/k) * (r(N+1) + h/2);
    
    % Lös systemet
    T_solve = A \ b;
    T_out(m) = T_solve(end);
    
    % Skapa hela temperaturprofilen (inklusive T0)
    T_full = [Ti; T_solve];
    
    % Plotta just denna profilen med rätt färg från färgskalan
    plot(r_full, T_full, 'Color', colors(m,:), 'LineWidth', 1.5);
end

% --- SNOFSA TILL DEN HÖGRA SUBPLOTTEN (Profilerna) ---
xlabel('Radie r (längdenheter)');
ylabel('Temperatur (^{\circ}C)');
title('Temperaturfördelning för olika \alpha');
grid on;
box on;
% Lägg till färgstapel (colorbar) för att visa vilket alpha varje linje har
colormap(parula);
c = colorbar;
c.Label.String = '\alpha (Värmeöverföringstal)';
clim([0 10]); % Används i nyare MATLAB istället för caxis

% --- RITA DEN VÄNSTRA SUBPLOTTEN (Uppgift 4d's huvudfråga) ---
subplot(1, 2, 1);
plot(alpha_vec, T_out, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
xlabel('Värmeöverföringstal \alpha');
ylabel('Temperatur vid ytterradien r=2 (^{\circ}C)');
title('Yttemperatur vs \alpha (N = 10000)');
grid on;
box on;