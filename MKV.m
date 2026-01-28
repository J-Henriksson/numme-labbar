clear
clc

graf_period = 10*365;
format long

% del a)

data = load("STHLMTEMP.mat");
T = data.Tdm(:,1);
t = (1:length(T))';

k = (2*pi)/365;
A = [ones(size(t)), sin(k*t), cos(k*t), sin(2*k*t), cos(2*k*t)];

c = A'*A \ A'*T

modell_1 = c(1) + c(2)*sin(k*t) + c(3)*cos(k*t) + ...
             c(4)*sin(2*k*t) + c(5)*cos(2*k*t);

figure(1)
clf
hold on
plot(modell_1(1:graf_period,1), "r")
plot(T(1:graf_period,1), "b")
title('mätdata och anpassad modell')
xlabel('tid (dagar)');
ylabel('Temperatur');
legend("modell", "mätdata")
hold off

% del b)

res_1 = T - modell_1;
MKsumma_1 = sum((res_1).^2)

figure(2)
clf
plot(res_1(1:graf_period,1))
title('Residual (T - modell_1)');
xlabel('tid (dagar)');
ylabel("differans")


%del c)

B = [ones(size(t)), t, t.^2, sin(k*t), cos(k*t), sin(2*k*t), cos(2*k*t)];
a = B'*B \ B'*T

modell_2 = a(1) + a(2)*t + a(3)*t.^2 + c(2)*sin(k*t) + ...
c(3)*cos(k*t) + c(4)*sin(2*k*t) + c(5)*cos(2*k*t);

figure(3)
clf
hold on
plot(modell_2(1:graf_period,1), "r")
plot(T(1:graf_period,1), "b")
title('mätdata och anpassad utökad modell')
xlabel('tid (dagar)');
ylabel('Temperatur');
legend("utökad modell", "mätdata")
hold off

% del d)

res_2 = T - modell_2;
MKsumma_2 = sum((res_2).^2)

figure(4)
clf
plot(res_2(1:graf_period,1))
title('Residual (T - modell_1)');
xlabel('tid (dagar)');
ylabel("differans")

% del e)

%medevärde första året
average_1 = sum(modell_2(1:365)) / 365;
%medelvärde andra året
average_2 = sum(modell_2(end - 364:end)) / 365;
%skillnad medelvärden
average_diff = average_2 - average_1

t_vand = -a(2) / (2*a(3));
ar_vanda = 1756 + (t_vand / 365.25);
fprintf('Trenden vänder ca år: %.0f\n', ar_vanda);

%{
analys av graf:
Genom att jämföra modellens medelvärde för det första och sista året ser vi att 
temperaturen har stigit med ca 1.4 grader totalt. 

Tittar vi på koefficienterna ser vi att ökningen inte har varit konstant. 
Eftersom a2-termen (framför t^2) är positiv så böjer kurvan sig uppåt, vilket 
betyder att uppvärmningen går snabbare och snabbare ju längre tiden går.

Det finns också en vändpunkt i modellen. Eftersom a1 är negativ börjar trenden 
med en svag sänkning, men eftersom a2 är positiv tar den kvadratiska termen 
över efter ett tag. Om man räknar på det (t = -a1/(2*a2)) ser man att det vänder 
runt år 1850. Efter det pekar trenden bara uppåt. Det här förklarar varför 
ökningen ser ut att ske "plötsligt" mot slutet av mätperioden jämfört med 
hur det såg ut på 1700-talet.
%}
 
