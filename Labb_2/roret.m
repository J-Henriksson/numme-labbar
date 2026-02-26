clear; clc; close all;

%% --- UPPGIFT 4c: Konvergensstudie för k = 1, alpha = 1 ---
fprintf('--- Uppgift 4c: Successiva fördubblingar ---\n');

% Fasta parametrar
Ti = 450;
Te = 20;
Ri = 1;
Ry = 2;
k = 1;
alpha_4c = 1;

N = 25;              % Startvärde enligt instruktion
T_out_old = 0;       % Variabel för att spara föregående resultat
skillnad = 1;        % Startvärde för loopen

% Loopa och dubbla N tills vi har stabiliserat minst 1 decimal (skillnad < 0.05)
while skillnad > 0.05
    h = (Ry - Ri) / (N + 1);
    r = Ri+h : h : Ry;
    
    A = sparse(N+1, N+1);
    b = zeros(N+1, 1);
    
    % Rad 1
    A(1, 1) = -2 * r(1); A(1, 2) = r(1) + h/2; b(1) = -(r(1) - h/2) * Ti;
    
    % Rad 2 till N
    for i = 2:N
        A(i, i-1) = r(i) - h/2;
        A(i, i)   = -2 * r(i);
        A(i, i+1) = r(i) + h/2;
    end
    
    % Sista raden
    A(N+1, N)   = 2 * r(N+1);
    A(N+1, N+1) = -2 * r(N+1) - (2*h*alpha_4c/k) * (r(N+1) + h/2);
    b(N+1)      = -(2*h*alpha_4c*Te/k) * (r(N+1) + h/2);
    
    % Lös systemet
    T_solve = A \ b;
    T_out_new = T_solve(end);
    
    % Beräkna skillnaden mot förra N för att kolla precisionen
    skillnad = abs(T_out_new - T_out_old);
    
    % Skriv ut i kommandofönstret
    fprintf('N = %4d ger temperaturen: %.5f\n', N, T_out_new);
    
    % Förbered för nästa varv om precisionen inte räcker
    if skillnad > 0.05
        T_out_old = T_out_new;
        N = N * 2; % Dubbla N
    end
end

fprintf('=> Kravet på en korrekt decimal uppnått vid N = %d\n\n', N);

% Plotta temperaturfördelningen för det sista (och bästa) N-värdet
figure(1);
T_full = [Ti; T_solve];
r_full = [Ri, r];
plot(r_full, T_full, 'b-', 'LineWidth', 1.5);
xlabel('Radie r (längdenheter)');
ylabel('Temperatur (^{\circ}C)');
title(['Uppgift 4c: Temperaturfördelning (N = ', num2str(N), ')']);
grid on;


%% --- UPPGIFT 4d: Variera alpha (N = 10000) ---

N_4d = 10000; % Nytt N enligt uppgift 4d
h = (Ry - Ri) / (N_4d + 1);
r = Ri+h : h : Ry;

% 20 intervall från 0 till 10 innebär 21 punkter
alpha_vec = linspace(0, 10, 21); 
T_out_4d = zeros(length(alpha_vec), 1);

% Bygg statiska delen av matrisen (Görs bara en gång för att spara tid!)
A = sparse(N_4d+1, N_4d+1);
b = zeros(N_4d+1, 1);

A(1, 1) = -2 * r(1); A(1, 2) = r(1) + h/2; b(1) = -(r(1) - h/2) * Ti;
for i = 2:N_4d
    A(i, i-1) = r(i) - h/2;
    A(i, i)   = -2 * r(i);
    A(i, i+1) = r(i) + h/2;
end

% Loopa över de olika alpha-värdena
for m = 1:length(alpha_vec)
    alpha_current = alpha_vec(m);
    
    % Byt bara ut sista raden
    A(N_4d+1, N_4d)   = 2 * r(N_4d+1);
    A(N_4d+1, N_4d+1) = -2 * r(N_4d+1) - (2*h*alpha_current/k) * (r(N_4d+1) + h/2);
    b(N_4d+1)         = -(2*h*alpha_current*Te/k) * (r(N_4d+1) + h/2);
    
    % Lös systemet och spara sista värdet
    T_solve = A \ b;
    T_out_4d(m) = T_solve(end);
end

% Plotta yttemperaturen som funktion av alpha
figure(2);
plot(alpha_vec, T_out_4d, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
xlabel('\alpha (Värmeöverföringstal)');
ylabel('Temperatur vid ytterradien r=2 (^{\circ}C)');
title('Uppgift 4d: Yttemperatur som funktion av \alpha (N = 10000)');
grid on;