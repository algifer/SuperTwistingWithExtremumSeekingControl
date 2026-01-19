%% Simulation of the LIPATA dark fermentation reactor for hidrogen production

clc; close all; clear all;

% File .xlsx: Latin Hypercube of data (Qin, GluIn, Hpr)
% File .mat: surface obtained from the data in the Excel file.
 

%% Particular case
filename_doe = 'IO_regression_By_DOE_with_Noise_5_noncumulated.xlsx';
filename_trained = 'trainedModel_with_5_sample_noncumulated.mat';


% Read data from excel file
data = readtable(filename_doe);

%%Load trained model
load(filename_trained);

% Get columns of GluIn and Qin
GluInLoad = data.GluIn;
QinLoad = data.Qin;

disp('GluIn:');
disp(GluInLoad);
disp('Qin:');
disp(QinLoad);


figure;
%% Subplot 1: Design of Experiments
subplot(2, 1, 1);
plot(GluInLoad, QinLoad, 'o', 'LineWidth', 2, 'MarkerSize', 6);
axis([0 30 0 6]);
title('a) Design of Experiment Qin vs. GluIn', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('GluIn [g/L]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Qin [L/d]', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 14, 'LineWidth', 1.5);
grid on;
box on;


% Variable at the reactor input
points = 80;
QinEq = linspace(1,5.5,points); %1 - 5 %2.5;            % L/d
GluInEq = linspace(5,30,points); %5 -30 %20.0;         % g/L

% Inputs Length
lQin = length(QinEq);
lGluIn = length(GluInEq);

Hpr_fit = zeros(lQin,lGluIn);

for j=1:lQin
    for k=1:lGluIn
        TI = table(GluInEq(k), QinEq(j),'VariableNames',{'GluIn', 'Qin'});
        
        Hpr_fit(j,k) = trainedModel.predictFcn(TI); % Particular case
    end
end


%% Subplot 2: Surface plot
subplot(2, 1, 2);
surf(GluInEq, QinEq, Hpr_fit);
xlabel('GluIn [g/L]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Qin [L/d]', 'FontSize', 14, 'FontWeight', 'bold');
zlabel('HPR [g(H_{2})/Ld]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');
title('b) Hydrogen Production Rate', 'FontSize', 14, 'FontWeight', 'bold');
colorbar;
set(gca, 'FontSize', 14, 'LineWidth', 1.5);
grid on;
box on;



