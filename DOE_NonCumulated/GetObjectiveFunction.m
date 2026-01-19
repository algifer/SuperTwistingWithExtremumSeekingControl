%% Simulation of the LIPATA dark fermentation reactor for hidrogen production

clc; close all; clear all;

% Especifica el nombre del archivo de Excel
filename_doe = 'IO_regression_By_DOE_without_Noise_35_main_points.xlsx';
%filename_doe = 'IO_regression_By_DOE_without_Noise_35.xlsx';
%filename_doe = 'IO_regression_By_DOE_without_Noise_35.xlsx';

filename_trained = 'trainedModel_without_35_sample_main_points.mat';
%filename_trained = 'trainedModel_without_35_sample_doe.mat';
%filename_trained = 'trainedModel_without_30_sample.mat';

% Lee los datos del archivo Excel en una tabla
data = readtable(filename_doe);

%%Load trained model
load(filename_trained);

%% Model in open loop
% Time definition
t0 = 0;   % d
tf = 120;    % d
Dt = 4;   % d

% Simulation time
t = t0:Dt/(24):tf;
lt = length(t);

% Inputs and Outputs initial values
Hpr_ode = zeros(1, lt);
Qin_ode = linspace(0,5.5,lt); %20.0;         % g/L
GluIn_ode = 25*ones(1, lt);

for k=1:lt-1
    % Equilibrium model
    TI = table(GluIn_ode(k), Qin_ode(k),'VariableNames',{'GluIn', 'Qin'});
    y = trainedModel.predictFcn(TI);
    Hpr_ode(k+1) = y; 
    
end
figure;
[maxHpr, pos] = max(Hpr_ode);
plot(t,Hpr_ode,'bl',t(pos),maxHpr,'*c','LineWidth',2)
%axis([0 30 0 35])
xlabel('time [d]');
ylabel('HPr [L(H_{2})L^{-1}d^{-1}]');

figure;
[maxHpr, pos] = max(Hpr_ode);
plot(Qin_ode,Hpr_ode,'bl',Qin_ode(pos),maxHpr,'*c','LineWidth',2)
axis([0 5.5 0 35])
xlabel('Input flow rate (Qin) [Ld^{-1}]');
ylabel('HPr [L(H_{2})L^{-1}d^{-1}]'); 