% Especifica el nombre del archivo de Excel
filename = 'IO_regression_By_DOE_with_Noise_35_main_points.xlsx';

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
plot(GluIn',Qin','o','LineWidth',2)
axis([0 35 0 7])
title('Desig of Experiment Qin vs GluIn')
xlabel('Input glucose (GluIn)[gL^{-1}]');
ylabel('Input flow rate (Qin) [Ld^{-1}]');
