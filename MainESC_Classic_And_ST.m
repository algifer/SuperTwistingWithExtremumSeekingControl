%Extremum Seeking Control


%% Simulation of the LIPATA dark fermentation reactor for hidrogen production

clc; close all; clear;

%% Load the model parameters and experimental data

% Mecanistic Model
parametersLIPATA2014;
dataLIPATA2014;


%% Time definition
t0 = 0;   % d
tf = 120;    % d
Dt = 4;   % h

% Simulation time
t = t0:Dt/(24):tf; %Sample time of 0.167 is each 4 hour
%t = t0:0.167:tf;
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
y0C = y0;

%% ODEs solution for Close Loop
%% Productivity in Close Loop ESC Super Twisting
%% Variable at the reactor input
HPR = zeros(1,lt);
Qin = 1.2*ones(1,lt);% Input flowrate (L/d)
Qin1 = zeros(1,lt); % Input flowrate (L/d)
% GluIn = 25*ones(1, lt);

GluIn = zeros(1, lt);
HPR_ref = zeros(1, lt);

for k=1:lt
    if((t(k)>=0) && (t(k)<5))
        GluIn(k) = 0;
        HPR_ref(k) = 0; 
    elseif((t(k)>=5) && (t(k)<25))
        GluIn(k) = 5;
        HPR_ref(k) = 5.3; 
    elseif((t(k)>=25) && (t(k)<50))
        GluIn(k) = 20;
        HPR_ref(k) = 25; 
    elseif((t(k)>=50) && (t(k)<75))
        GluIn(k) = 25;
        HPR_ref(k) = 32; 
    elseif((t(k)>=75) && (t(k)<95))
        GluIn(k) = 10;
        HPR_ref(k) = 11.4;
    else
        GluIn(k) = 15;
        HPR_ref(k) = 18.45;
    end
end

% figure;
% plot(t, HPR_ref, '--r', 'LineWidth',2)
% xlabel('Time (d)');
% ylabel('HPR-REF (L[H_2]L^{-1}d^{-1})');


y_ = zeros(10,lt);
y_(:,1) = y0;

HPR(1) = 0;
gradEst = ones(1,lt);

% Optimization beginning
tic = 5;                % d

% Extremum seeking parameters
% alphaST = 0.37; %0.4; %0.0725;
% lambdaST = 0.1; %0.091; %0.1;

alphaST = 0.4; %0.0725;
lambdaST = 0.091; %0.1;

% Integrator options
optionsST = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
optionsInt = odeset('RelTol',1e-6,'AbsTol',1e-6);
for p=1:lt-1    
    % Trained Model with Noise
    [tOut3, yOut] = ode15s('H2ReactorLIPATA2RODE', [t(p) t(p+1)], y0, optionsST, GluIn(p), Qin(p), V, Vgas, K, ...
                      H2MM, muMax, KGlu, KHCO2, KHH2, KlaCO2, KlaH2, Patm, pVapH2O, R, TReac, TAmb);

    y_(:,p+1) = yOut(end,:)';
    y0 = y_(:,p+1);

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
    y = QH2gas/V;
    HPR(p+1) = y; 
    
    if t(p)>=tic
        fprintf('\n Productivity=%f, Inflow=%f, Gradient=%f',HPR(p),Qin(p-1),gradEst(p));
        
        % Gradient approximation
        if t(p)>tic
            gradEst(p+1) = (HPR(p+1)-HPR(p))/(Qin(p)-Qin(p-1));
        end
        
        % Optimization algorithm
        [~, QinOut] = ode15s(@integratorST, [t(p) t(p+1)], Qin1(p), optionsInt,...
                            alphaST, gradEst(p));
        Qin1(p+1) = QinOut(end);
        
        Qin(p+1) = (lambdaST*sqrt(abs(gradEst(p+1)))*sign(gradEst(p+1))) + Qin1(p+1);
    end
end


%% ODEs solution for Close Loop
%% Productivity in Close Loop ESC Classic
%% Variable at the reactor input
HPR_C = zeros(1, lt);
Qin_C = 1.2*ones(1, lt);
% GluIn_C = 25*ones(1, lt);
GluIn_C = GluIn;

y_C = zeros(10,lt);
y_C(:,1) = y0C;

HPR_C(1) = 0;
gradEst_C = 1*ones(1,lt);

% Integrator options
optionsC = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

%% Initial conditions
Eta0 = 0;
Qin_est0 = 0;

y1C = Eta0;
y2C = Qin_est0;

w = 120; %rad/day
A = 0.3;
wh = 5;
kI = 4;


    for k=1:lt-1

        % Call to the IVP solver
        [tOutC, yOutC] = ode15s('H2ReactorLIPATA2RODE', [t(k) t(k+1)], y0C, optionsC, GluIn_C(k), Qin_C(k), V, Vgas, K, ...
                          H2MM, muMax, KGlu, KHCO2, KHH2, KlaCO2, KlaH2, Patm, pVapH2O, R, TReac, TAmb);
    
        y_C(:,k+1) = yOutC(end,:)';
        y0C = y_C(:,k+1);
        %% Plot results
        Glu_C = yOutC(end,1);           % g/L
        Ace_C = yOutC(end,2);           % g/L
        Pro_C = yOutC(end,3);           % g/L
        Bu_C = yOutC(end,4);            % g/L
        EtOH_C = yOutC(end,5);          % g/L
        X_C = yOutC(end,6);             % g/L
        CO2_C = yOutC(end,7);           % M/L
        H2_C = yOutC(end,8);            % g/L
        CO2gas_C = yOutC(end,9);        % M/L
        H2gas_C = yOutC(end,10);        % g/L

        pGasH2_C = (H2gas_C*R*TReac)/H2MM;
        rhoTH2_C = KlaH2 * (H2_C-(H2MM*KHH2*pGasH2_C));
        QH2gas_C = ((R*TAmb)/(Patm-pVapH2O))*V*(rhoTH2_C/H2MM);
        yC = QH2gas_C/V;
        HPR_C(k+1) = yC;
        
        if t(k)>=tic
            fprintf('\n Productivity=%f, Inflow=%f, Gradient=%f',HPR_C(k),Qin_C(k-1),gradEst_C(k));
        
            % Gradient approximation
            [tOut2C, eta] = ode15s(@HighPassFilterESC_ODE, [t(k) t(k+1)], y1C, optionsC, yC, wh);
            y1C = eta(end,1);
            cita = (yC - y1C)*A*sin(w*t(k+1));
            gradEst_C(k+1) = cita;
            
            
            [tOut3C, Qin_est] = ode15s(@IntegratorESC_ODE, [t(k) t(k+1)], y2C, optionsC, cita, kI);
            y2C = Qin_est(end,1);
            Qin_C(k+1) = y2C + A*sin(w*t(k+1));
        end

        
        
    end


%% Plot results
% figure;
% plot(t, Qin_C, t, Qin, 'LineWidth',2);
% xlabel('Time (d)');
% ylabel('Qin (Ld^{-1})');
% %title('Input flow rate')
% legend('Qin with Classic ESC', 'Qin with Super Twisting ESC')
% 
% figure;
% plot(t, HPR_C, t, HPR, 'LineWidth',2)
% xlabel('Time (d)');
% ylabel('HPR (L[H_2]L^{-1}d^{-1})');
% %title('Hydrogen Production Rate')
% legend('HPR with Classic ESC', 'HPR with Super Twisting ESC')
% 
% figure;
% plot(t, gradEst_C, t, gradEst, 'LineWidth',2);
% xlabel('Time (d)');
% ylabel('\nabla HPR');
% %title('Hydrogen Production Rate Gradient')
% legend('\nabla HPR with Classic ESC', '\nabla HPR with Super Twisting ESC')


%% Plot results
figure;
plot(t, GluIn, 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Time [d]', 'FontSize', 12, 'FontWeight', 'bold');
%ylabel('Glu_{in} (Ld^{-1})');
ylabel('GluIn [g/d]', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0 30])
%title('Input flow rate')
%legend('GluIn with Super Twisting ESC')
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;
box on;


figure;
plot(t, Qin, 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Time [d]', 'FontSize', 12, 'FontWeight', 'bold');
%ylabel('Qin (Ld^{-1})');
ylabel('Qin [L/d]', 'FontSize', 12, 'FontWeight', 'bold');
%title('Input flow rate')
%legend('Qin with Super Twisting ESC')
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;
box on;

figure;
plot(t, HPR_ref, '--r', 'LineWidth', 2, 'MarkerSize', 6)
hold on
plot(t, HPR, 'LineWidth', 2, 'MarkerSize', 6)
xlabel('Time [d]', 'FontSize', 12, 'FontWeight', 'bold');
%ylabel('HPR (L[H_2]L^{-1}d^{-1})');
ylabel('HPR [g[H_2]/Ld]', 'FontSize', 12, 'FontWeight', 'bold');
%title('Hydrogen Production Rate')
legend('Off-line Maximum HPR', 'HPR Super Twisting based ESC', 'FontSize', 10)
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
box on;
hold off

figure;
plot(t, gradEst, 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Time [d]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\nabla HPR', 'FontSize', 12, 'FontWeight', 'bold');
%title('Hydrogen Production Rate Gradient')
%legend('\nabla HPR with Super Twisting ESC')
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
grid on;
box on;