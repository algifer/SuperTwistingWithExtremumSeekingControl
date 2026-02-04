%Extremum Seeking Control


%% Simulation of the LIPATA dark fermentation reactor for hidrogen production

clc; close all; clear;

%% Load the model parameters and experimental data

% Mecanistic Model
parametersLIPATA2014;
dataLIPATA2014;

% Trained Model With Noise
load('trainedModel_with_15_sample_cumulated.mat');


%% Time definition
t0 = 0;   % d
tf = 60;    % d
Dt = 4;   % h u


% Simulation time
t = t0:Dt/(24):tf; %Sample time of 0.167 is each 4 hour
%t = t0:0.167:tf;
lt = length(t);


%% Variable at the reactor input
%% Response surface
points = 80;
Qin_Open_Loop = linspace(1,5.5,points); %1 - 5 %2.5;            % L/d
GluIn_Open_Loop = linspace(5,30,points); %5 -30 %20.0;         % g/L

% Inputs Length
lQin = length(Qin_Open_Loop);
lGluIn = length(GluIn_Open_Loop);

Hpr_fit_Open_Loop = zeros(lQin,lGluIn);

for j=1:lQin
    for k=1:lGluIn
        TI = table(GluIn_Open_Loop(k), Qin_Open_Loop(j),'VariableNames',{'GluIn', 'Qin'});
        %Hpr_fit_Open_Loop(j,k) = trainedModelWithNoise.predictFcn(TI);
        Hpr_fit_Open_Loop(j,k) = trainedModel.predictFcn(TI);
    end
end

figure;
surf(GluIn_Open_Loop,Qin_Open_Loop,Hpr_fit_Open_Loop)
ylabel('Q_{in} (L/d)');
xlabel('Glu_{in} (g/L)');
zlabel('HPR (g[H_{2}]/Ld)');
%title('HPR Response Surface')
colorbar

%% Productivity in open Loop
Hpr_open_loop = zeros(1,lt);
Qin_open_loop = linspace(0,5.34,lt); %20.0;         % g/L
GluIn_open_loop = 25*ones(1,lt);

for k=1:lt-1
    % Trained Model with Noise
    TI = table(GluIn_open_loop(k), Qin_open_loop(k),'VariableNames',{'GluIn', 'Qin'});
    %y_ol = trainedModelWithNoise.predictFcn(TI);
    y_ol = trainedModel.predictFcn(TI);
    Hpr_open_loop(k+1) = y_ol;
end


figure;
[maxHpr_ol, pos] = max(Hpr_open_loop);
plot(t,Hpr_open_loop,t(pos),maxHpr_ol,'*c','LineWidth',2)
%title('Hydrogen production rate in Open Loop with GluIn = 25')
xlabel('Time (d)');
% ylabel('HPR (L(H_{2})L^{-1}d^{-1})');
ylabel('HPR (g[H_{2}]/Ld)');

%% Initial conditions
Glu0 = GluVec(1);            % g/L
Ace0 = AceVec(1);            % g/L
Pro0 = ProVec(1);            % g/L
Bu0 = BuVec(1);              % g/L
EtOH0 = EtOHVec(1);          % g/L
X0 = XVec(1);                % g/L

qCO2Vgas0 = qCO2VgasVec(1);     % L/d
qH2Vgas0 = qH2VgasVec(1);       % L/d

qCO2gas0 = ((Patm-pVapH2O)*qCO2Vgas0)/(R*TAmb*V);                   % mol/Ld
qH2gas0 = ((Patm-pVapH2O)*qH2Vgas0)/(R*TAmb*V);                     % mol/Ld

pCO2gas0 = (qCO2Vgas0/(qCO2Vgas0+qH2Vgas0)) * (Patm-pVapH2O);       % bar
CO2gas0 = pCO2gas0/(R*TReac);                                       % mol/L
CO20 = (qCO2gas0/KlaCO2) + (KHCO2*pCO2gas0);                        % mol/L

pH2gas0 = (qH2Vgas0/(qCO2Vgas0+qH2Vgas0)) * (Patm-pVapH2O);         % bar
H2gas0 = (pH2gas0/(R*TReac))*H2MM;                                  % g/L
H20 = (qH2gas0/KlaH2) + (KHH2*pH2gas0*H2MM);                        % g/L

y0 = [Glu0; Ace0; Pro0; Bu0; EtOH0; X0; CO20; H20; CO2gas0; H2gas0];

y_ = zeros(10,lt);
y_(:,1) = y0;

%% ODEs solution
HPR_mechanistic = zeros(1,lt);
HPR_interpolated = zeros(1,lt);
HPR_interpolated_2 = zeros(1,lt);
Qin_mechanistic = zeros(1,lt);
Qin1_mechanistic = zeros(1,lt);
Qin_interpolated = zeros(1,lt);
Qin1_interpolated = zeros(1,lt);
Qin_interpolated_2 = zeros(1,lt);
Qin1_interpolated_2 = zeros(1,lt);
GluIn = 25*ones(1, lt);

HPR_mechanistic(1) = 0;
HPR_interpolated(1) = 0;
HPR_interpolated_2(1) = 0;

gradEst_mechanistic = ones(1,lt);
gradEst_interpolated = ones(1,lt);
gradEst_interpolated_2 = ones(1,lt);

% Optimization beginning
tic = 5;                % d

% Extremum seeking parameters
% alphaST = 0.37; %0.4; %0.0725;
% lambdaST = 0.1; %0.091; %0.1;

flag        = 1;
alphaST = 0.4;
lambdaST = 0.091;

options_mechanistic = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
optionsInt_mechanistic = odeset('RelTol',1e-6,'AbsTol',1e-6);
optionsInt_interpolated = odeset('RelTol',1e-6,'AbsTol',1e-6);
optionsInt_interpolated_2 = odeset('RelTol',1e-6,'AbsTol',1e-6);

for k=1:lt-1 
    
    % Mechanistic Model
    [tOut3, yOut] = ode15s('H2ReactorLIPATA2RODE', [t(k) t(k+1)], y0, options_mechanistic, GluIn(k), Qin_mechanistic(k), V, Vgas, K, ...
                      H2MM, muMax, KGlu, KHCO2, KHH2, KlaCO2, KlaH2, Patm, pVapH2O, R, TReac, TAmb);

    y_(:,k+1) = yOut(end,:)';
    y0 = y_(:,k+1);

    %% Plot results
    Glu = yOut(end,1);           % g/L
    Ace = yOut(end,2);           % g/L
    Pro = yOut(end,3);           % g/L
    Bu = yOut(end,4);            % g/L
    EtOH = yOut(end,5);          % g/L
    X = yOut(end,6);             % g/L
    CO2 = yOut(end,7);           % M/L
    H2 = yOut(end,8);            % g/L
    CO2gas = yOut(end,9);        % M/L
    H2gas = yOut(end,10);        % g/L

    pGasH2 = (H2gas*R*TReac)/H2MM;
    rhoTH2 = KlaH2 * (H2-(H2MM*KHH2*pGasH2));
    QH2gas = ((R*TAmb)/(Patm-pVapH2O))*V*(rhoTH2/H2MM);
    y1 = QH2gas/V;
    
    
    if t(k)>=tic
        HPR_mechanistic(k+1) = y1; 
        fprintf('\n Mechanistic HPR=%f, Inflow=%f, Gradient=%f',HPR_mechanistic(k),Qin_mechanistic(k-1),gradEst_mechanistic(k));
        
        % Gradient approximation
        if t(k)>tic
            gradEst_mechanistic(k+1) = (HPR_mechanistic(k+1)-HPR_mechanistic(k))/(Qin_mechanistic(k)-Qin_mechanistic(k-1));
        end
        
        % Optimization algorithm
        [~, QinOut1] = ode15s(@integratorST, [t(k) t(k+1)], Qin1_mechanistic(k), optionsInt_mechanistic,...
                            alphaST, gradEst_mechanistic(k));
        Qin1_mechanistic(k+1) = QinOut1(end);
        
        Qin_mechanistic(k+1) = (lambdaST*sqrt(abs(gradEst_mechanistic(k+1)))*sign(gradEst_mechanistic(k+1))) + Qin1_mechanistic(k+1);
    end
    
    

%     %% No Acumulated Trained Model
%     TI = table(GluIn(k), Qin_interpolated(k),'VariableNames',{'GluIn', 'Qin'});
%     y2 = trainedModelWithNoise.predictFcn(TI);
%     
%     if t(k)>=tic
%         HPR_interpolated(k+1) = y2; 
%         fprintf('\n Interpolated HPR=%f, Inflow=%f, Gradient=%f',HPR_interpolated(k),Qin_interpolated(k-1),gradEst_interpolated(k));
%         
%         % Gradient approximation
%         if t(k)>tic
%             gradEst_interpolated(k+1) = (HPR_interpolated(k+1)-HPR_interpolated(k))/(Qin_interpolated(k)-Qin_interpolated(k-1));
%         end
%         
%         % Optimization algorithm
%         [~, QinOut2] = ode15s(@integratorST, [t(k) t(k+1)], Qin1_interpolated(k), optionsInt_interpolated,...
%                             alphaST, gradEst_interpolated(k));
%         Qin1_interpolated(k+1) = QinOut2(end);
%         
%         Qin_interpolated(k+1) = (lambdaST*sqrt(abs(gradEst_interpolated(k+1)))*sign(gradEst_interpolated(k+1))) + Qin1_interpolated(k+1);
%     end
%     
    
    %% Acumulated Trained Model
    TI2 = table(GluIn(k), Qin_interpolated_2(k),'VariableNames',{'GluIn', 'Qin'});
    y3 = trainedModel.predictFcn(TI2);
    
    if t(k)>=tic
        HPR_interpolated_2(k+1) = y3; 
        fprintf('\n Interpolated HPR=%f, Inflow=%f, Gradient=%f\n',HPR_interpolated_2(k),Qin_interpolated_2(k-1),gradEst_interpolated_2(k));
        
        % Gradient approximation
        if t(k)>tic
            gradEst_interpolated_2(k+1) = (HPR_interpolated_2(k+1)-HPR_interpolated_2(k))/(Qin_interpolated_2(k)-Qin_interpolated_2(k-1));
        end
        
        % Optimization algorithm
        [~, QinOut_2] = ode15s(@integratorST, [t(k) t(k+1)], Qin1_interpolated_2(k), optionsInt_interpolated_2,...
                            alphaST, gradEst_interpolated_2(k));
        Qin1_interpolated_2(k+1) = QinOut_2(end);
        
        Qin_interpolated_2(k+1) = (lambdaST*sqrt(abs(gradEst_interpolated_2(k+1)))*sign(gradEst_interpolated_2(k+1))) + Qin1_interpolated_2(k+1);
    end
end

% %% Plot results
% figure;
% plot(t, Qin_mechanistic, 'LineWidth', 2, 'MarkerSize', 6);
% hold on
% plot(t, Qin_interpolated, t, Qin_interpolated_2, 'LineWidth', 2, 'MarkerSize', 6);
% xlabel('Time [d]', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Qin [L/d]', 'FontSize', 12, 'FontWeight', 'bold');
% %title('Input flow rate')
% legend('Qin to Mechanistic Model', ...
%        'Qin to Model without Cumulated Points', ...
%        'Qin to Model with Cumulated Points', 'FontSize', 10)
% set(gca, 'FontSize', 12, 'LineWidth', 1.5)
% grid on
% box on
% 
% figure;
% plot(t, HPR_mechanistic, 'LineWidth', 2, 'MarkerSize', 6)
% hold on
% plot(t, HPR_interpolated, t, HPR_interpolated_2, 'LineWidth', 2, 'MarkerSize', 6)
% xlabel('Time [d]', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('HPR [g[H_2]/Ld]', 'FontSize', 12, 'FontWeight', 'bold');
% %title('Hydrogen Production Rate')
% legend('HPR to Mechanistic Model', ...
%        'HPR to Model without Cumulated Points', ...
%        'HPR to Model with Cumulated Points', 'FontSize', 10)
% set(gca, 'FontSize', 12, 'LineWidth', 1.5)
% grid on
% box on
% 
% figure;
% plot(t, gradEst_mechanistic, 'LineWidth', 2, 'MarkerSize', 6);
% hold on
% plot(t, gradEst_interpolated, t, gradEst_interpolated_2, 'LineWidth', 2, 'MarkerSize', 6);
% xlabel('Time [d]', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('\nabla HPR', 'FontSize', 12, 'FontWeight', 'bold');
% %title('Hydrogen Production Rate Gradient')
% legend('Gradient Mechanistic Model', ...
%        'Gradient Model without Cumulated Points', ...
%        'Gradient Model with Cumulated Points', 'FontSize', 10)
% set(gca, 'FontSize', 12, 'LineWidth', 1.5)
% grid on
% box on


%% ============================================================
%% Numerical Hessian analysis for concavity and unimodality
%% (Reviewer requested quantitative verification)
%% ============================================================

% Same operating domain used in the response surface analysis
Qin_vec   = linspace(1,5.5,80);     % L/d
GluIn_vec = linspace(5,30,80);      % g/L

% Hessiano en los puntos optimos a lo largo de la superficie (negativo)
% Gradiente debe ser cero en los puntos optimos
% ERMS entre la diferencia de los puntos interpolados y los puntos muestreados

dQ = Qin_vec(2)   - Qin_vec(1);
dG = GluIn_vec(2) - GluIn_vec(1);

[Qgrid, Ggrid] = meshgrid(Qin_vec, GluIn_vec);

HPR_map = zeros(size(Qgrid));

% Evaluate static HPR map using the trained model
for i = 1:numel(Qgrid)
    TI = table(Ggrid(i), Qgrid(i), ...
        'VariableNames', {'GluIn','Qin'});
    %HPR_map(i) = trainedModelWithNoise.predictFcn(TI);
    HPR_map(i) = trainedModel.predictFcn(TI);
end

%% ============================================================
%% Optimal Qin for each GluIn and local curvature analysis
%% ============================================================

nG = length(GluIn_vec);

Qin_opt_vec   = zeros(nG,1);
HPR_opt_vec   = zeros(nG,1);
grad_Q_opt    = zeros(nG,1);
hess_Q_opt    = zeros(nG,1);
rmse_quad     = zeros(nG,1);

for iG = 1:nG

    % --- HPR profile for fixed GluIn
    HPR_profile = HPR_map(iG, :);

    % --- Locate optimal Qin (discrete maximum)
    [HPR_opt, iQ] = max(HPR_profile);

    Qin_opt_vec(iG) = Qin_vec(iQ);
    HPR_opt_vec(iG) = HPR_opt;

    % --- Gradient and Hessian w.r.t Qin
    if iQ > 1 && iQ < length(Qin_vec)

        % Central differences
        grad_Q_opt(iG) = (HPR_profile(iQ+1) - HPR_profile(iQ-1)) / (2*dQ);
        hess_Q_opt(iG) = (HPR_profile(iQ+1) ...
                         - 2*HPR_profile(iQ) ...
                         + HPR_profile(iQ-1)) / (dQ^2);

        % --- Local quadratic approximation error (RMSE)
        Qin_local = Qin_vec(iQ-1:iQ+1);
        HPR_local = HPR_profile(iQ-1:iQ+1);

        p = polyfit(Qin_local, HPR_local, 2);
        HPR_quad = polyval(p, Qin_local);

        rmse_quad(iG) = sqrt(mean((HPR_local - HPR_quad).^2));

    else
        % Boundary optima (skip curvature test)
        grad_Q_opt(iG) = NaN;
        hess_Q_opt(iG) = NaN;
        rmse_quad(iG)  = NaN;
    end
end

%% ============================================================
%% Figures requested by the reviewer
%% ============================================================

% --- GluIn vs Gradient
figure;
plot(GluIn_vec, grad_Q_opt, 'LineWidth', 2)
xlabel('Glu_{In} (g/L)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('\partial HPR / \partial Q_{in}', 'FontSize', 12, 'FontWeight', 'bold')
title('Gradient at optimal Q_{in}')
grid on
box on

% --- GluIn vs Hessian
figure;
plot(GluIn_vec, hess_Q_opt, 'LineWidth', 2)
xlabel('Glu_{In} (g/L)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('\partial^2 HPR / \partial Q_{in}^2', 'FontSize', 12, 'FontWeight', 'bold')
title('Local Hessian at optimal Q_{in}')
grid on
box on



% Óptimo global real de la superficie
[HPR_global_opt, idx] = max(HPR_map(:));
[iGopt, iQopt] = ind2sub(size(HPR_map), idx);

GluIn_opt = GluIn_vec(iGopt);
Qin_opt   = Qin_vec(iQopt);

% Paso pequeño coherente con la malla
dQ = Qin_vec(2) - Qin_vec(1);

% Gradiente REAL del modelo entrenado
HPR_plus  = trainedModel.predictFcn( ...
    table(GluIn_opt, Qin_opt + dQ, ...
    'VariableNames', {'GluIn','Qin'}) );

HPR_minus = trainedModel.predictFcn( ...
    table(GluIn_opt, Qin_opt - dQ, ...
    'VariableNames', {'GluIn','Qin'}) );

grad_real = (HPR_plus - HPR_minus) / (2*dQ);
grad_norm = abs(grad_real);

fprintf('||?HPR|| = %.4e\n',grad_norm);


