clc, close all, clear all;

%% Time definition
t0 = 0;   % d
tf = 60;    % d
Dt = 4;   % h

% Simulation time
t = t0:Dt/(24):tf; %Sample time of 0.167 is each 4 hour
%t = t0:0.167:tf;
lt = length(t);

GluIn = zeros(1, lt);  % Inicializar el vector GluIn con ceros

% Definir los puntos de cambio en el tiempo (puedes ajustarlos según lt)
changePoints = [1, floor(lt/5), floor(2*lt/5), floor(3*lt/5), floor(4*lt/5), lt];

% Asignar los valores correspondientes en cada intervalo
GluIn(changePoints(1):changePoints(2)) = 10;
GluIn(changePoints(2)+1:changePoints(3)) = 20;
GluIn(changePoints(3)+1:changePoints(4)) = 5;
GluIn(changePoints(4)+1:changePoints(5)) = 30;
GluIn(changePoints(5)+1:changePoints(6)) = 25;

%% Plot results
figure;
plot(t, GluIn, 'LineWidth',2);
xlabel('Time (d)');
ylabel('GluIn (Ld^{-1})');
%title('Input flow rate')
legend('GluIn with Super Twisting ESC')


