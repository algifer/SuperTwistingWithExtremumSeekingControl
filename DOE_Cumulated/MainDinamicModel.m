
%% Simulation of the LIPATA dark fermentation reactor for hidrogen production

clc; close all; clear;

%% Load the model parameters and experimental data

% Mecanistic Model
parametersLIPATA2014;
dataLIPATA2014;


%% Time definition
t0 = 0;   % d
tf = 120;    %60% d
Dt = 4;   % h

% Simulation time
t = t0:Dt/(24):tf; %Sample time of 0.167 is each 4 hour
lt = length(t);


%% Variable at the reactor input
%% Initial conditions;
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


%% ODEs solution response surface
%% Productivity in open Loop
%% Variable at the reactor input
%Qin_open_loop = linspace(0,5.34,lt); %20.0;         % g/L
Qin_open_loop = 4.2*ones(1,lt);                        % g/L
GluIn_open_loop = 20.83*ones(1,lt);


Hpr_open_loop = zeros(1,lt);

Xout = zeros(10,lt);
y_ = zeros(10,lt);
y_(:,1) = y0;

% Integrator options
options1 = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
steps = lt-1;


for p=1:steps

    % Call to the IVP solver
    [tOut2, yOut2] = ode15s('H2ReactorLIPATA2RODE', t, y0, options1, GluIn_open_loop(p), Qin_open_loop(p), V, Vgas, K, ...
                          H2MM, muMax, KGlu, KHCO2, KHH2, KlaCO2, KlaH2, Patm, pVapH2O, R, TReac, TAmb);

%     %% Plot results
%     Glu = yOut2(:,1);           % g/L
%     X = yOut2(:,2);             % g/L
%     H2 = yOut2(:,3);            % g/L
%     H2gas = yOut2(:,4);        % g/L
% 
%     pGasH2 = (H2gas*R*TReac)/H2MM;
%     rhoTH2 = KlaH2 * (H2-(H2MM*KHH2*pGasH2));
%     QH2gas = ((R*TAmb)/(Patm-pVapH2O))*V*(rhoTH2/H2MM);
%     %Hpr_open_loop(p) = QH2gas./V;
%     Hpr_open_loop = QH2gas./V;
%     
% %     Xout(1,p) = Glu;
% %     Xout(2,p) = X;
% %     Xout(3,p) = H2;
% %     Xout(4,p) = H2gas;

    y_(:,p+1) = yOut2(end,:)';
    y0 = y_(:,p+1);

    %% Plot results
    Glu = yOut2(end,1);           % g/L
    Ace = yOut2(end,2);           % g/L
    Pro = yOut2(end,3);           % g/L
    Bu = yOut2(end,4);            % g/L
    EtOH = yOut2(end,5);          % g/L
    X = yOut2(end,6);             % g/L
    CO2 = yOut2(end,7);           % M/L
    H2 = yOut2(end,8);            % g/L
    CO2gas = yOut2(end,9);        % M/L
    H2gas = yOut2(end,10);        % g/L

    pGasH2 = (H2gas*R*TReac)/H2MM;
    rhoTH2 = KlaH2 * (H2-(H2MM*KHH2*pGasH2));
    QH2gas = ((R*TAmb)/(Patm-pVapH2O))*V*(rhoTH2/H2MM);
    y = QH2gas./V;
    Hpr_open_loop = y; 

    Xs = ['Qin: ', num2str(Qin_open_loop(p)),' GluIn: ',num2str(GluIn_open_loop(p)),' H: ', num2str(Hpr_open_loop(p))];
    disp(Xs)
    break

end

%% Plot results
%Inputs
% figure(1), plot(t, Qin_open_loop,'LineWidth',2), xlabel('Time (d)'), ylabel('Input flow rate (Ld^{-1})');
% figure(2), plot(t, GluIn_open_loop,'LineWidth',2), xlabel('Time (d)'), ylabel('Input Glu  (Ld^{-1})');

%States
% figure(3), plot(t,Glu,'LineWidth',2), xlabel('Time (d)'), ylabel('Glucose (g/L)');
% figure(4), plot(t,X,'LineWidth',2), xlabel('Time (d)'), ylabel('Biomass (g/L)');
% figure(5), plot(t,H2,'LineWidth',2), xlabel('Time (d)'), ylabel('H2');
% figure(6), plot(t,H2gas,'LineWidth',2), xlabel('Time (d)'), ylabel('H2 gas flow (mL/h)');

% % figure(3), plot(t,Xout(1,:)), xlabel('Time (d)'), ylabel('Glucose (g/L)');
% % figure(4), plot(t,Xout(2,:),'-b'), xlabel('Time (d)'), ylabel('Biomass (g/L)');
% % figure(5), plot(t,Xout(3,:)), xlabel('Time (d)'), ylabel('H2');
% % figure(6), plot(t,Xout(4,:)), xlabel('Time (d)'), ylabel('H2 gas flow (mL/h)');

% Output
figure(7), plot(t, Hpr_open_loop,'LineWidth',2), xlabel('Time (d)'), ylabel('Hydrogen production rate (L[H_2]L^{-1}d^{-1})');

