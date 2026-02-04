clc; close all; clear all;
% Especifica el nombre del archivo de Excel
%filename = 'IO_regression_By_DOE_with_Noise_35_main_points.xlsx';
filename = 'IO_regression_By_DOE_with_Noise_30_noncumulated.xlsx';

% Lee los datos del archivo Excel en una tabla
data = readtable(filename);

% Extrae las columnas GluIn, Qin, y Hpr
GluIn = data.GluIn;
Qin = data.Qin;
Hpr = data.Hpr;

% Mostrar los datos en la consola (opcional)
disp('GluIn:');
disp(GluIn);
disp('Qin:');
disp(Qin);
disp('Hpr:');
disp(Hpr);

% Graficar los resultados
figure(1);
%plot(GluIn',Qin','o','LineWidth',2)
plot(GluIn',Qin','o','LineWidth',1)
%plot(GluIn',Qin','o')
axis([0 35 0 7])
title('Desig of Experiment Qin vs GluIn')
xlabel('Input glucose (GluIn)[gL^{-1}]');
ylabel('Input flow rate (Qin) [Ld^{-1}]');

% Gráfica 3D
figure(2);
scatter3(GluIn, Qin, Hpr, 60, Hpr, 'filled','LineWidth',1)
grid on
box on

xlabel('w','FontSize', 14, 'FontWeight', 'bold');
ylabel('u','FontSize', 14, 'FontWeight', 'bold');
zlabel('y','FontSize', 14, 'FontWeight', 'bold');

% xlabel('GluIn [g/L]');
% ylabel('Qin [L/d]');
% zlabel('HPR [g(H_{2})/Ld]');
%title('Data')

% colorbar
% view(135,30)
