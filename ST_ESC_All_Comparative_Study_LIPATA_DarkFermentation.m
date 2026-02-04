%Extremum Seeking Control


%% Simulation of the LIPATA dark fermentation reactor for hidrogen production

clc; close all; clear;

%% Load the model parameters and experimental data

% Mecanistic Model
parametersLIPATA2014;
dataLIPATA2014;

% Trained Model With Noise
load('trainedModel_with_30_sample_best_selected_noncumulated.mat'); % No acumulated
load('trainedModel_with_30_sample_cumulated.mat'); % acumulated


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
        Hpr_fit_Open_Loop(j,k) = trainedModelWithNoise.predictFcn(TI);
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
    y_ol = trainedModelWithNoise.predictFcn(TI);
    Hpr_open_loop(k+1) = y_ol;
end


% figure;
% [maxHpr_ol, pos] = max(Hpr_open_loop);
% plot(t,Hpr_open_loop,t(pos),maxHpr_ol,'*c','LineWidth',2)
% %title('Hydrogen production rate in Open Loop with GluIn = 25')
% xlabel('Time (d)');
% % ylabel('HPR (L(H_{2})L^{-1}d^{-1})');
% ylabel('HPR (g[H_{2}]/Ld)');

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
Qin_mechanistic = zeros(1,lt);
Qin1_mechanistic = zeros(1,lt);
Qin_interpolated = zeros(1,lt);
Qin1_interpolated = zeros(1,lt);

HPR_mechanistic(1) = 0;
HPR_interpolated(1) = 0;
HPR_interpolated_2(1) = 0;

gradEst_mechanistic = ones(1,lt);
gradEst_interpolated = ones(1,lt);
gradEst_interpolated_2 = ones(1,lt);

options_mechanistic = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
optionsInt_mechanistic = odeset('RelTol',1e-6,'AbsTol',1e-6);
optionsInt_interpolated = odeset('RelTol',1e-6,'AbsTol',1e-6);
optionsInt_interpolated_2 = odeset('RelTol',1e-6,'AbsTol',1e-6);

%% Storage variables
HPR_PI  = zeros(1,lt);
HPR_ESC  = zeros(1,lt);
HPR_STESC = zeros(1,lt);

% Qin_PID  = 1.5*ones(1,lt);
% Qin_ESC  = 1.5*ones(1,lt);
% Qin_STESC = 1.5*ones(1,lt);

Qhat_ESC = zeros(1, lt);

Qin_PI  = 0.1*zeros(1,lt);
Qin_ESC  = 0.1*zeros(1,lt);
Qin_STESC = 0.1*zeros(1,lt);

GluIn = 25*ones(1,lt);

%% Controller parameters
tic = 5;                % d
flag = 1;

% PI
Kp = 0.15; Ki = 0.025; e_int = 0;

% Classical ESC
w = 120; %rad/day
A = 0.3;
wh = 5;
kI = 4;

Qmin = 0.5;
Qmax = 5.5;

y2C = 0;
y1C = 0;

% ST-ESC (proposed)
alphaST = 0.4; lambdaST = 0.091;
gradST = ones(1,lt); Qin1 = zeros(1,lt);

options = odeset('RelTol',1e-6,'AbsTol',1e-6);


%% Golden search (offline reference)
Qtest = linspace(0.5,5.5,80);
HPR_static = zeros(size(Qtest));
for i=1:length(Qtest)
    TI = table(25,Qtest(i),'VariableNames',{'GluIn','Qin'});
    HPR_static(i) = trainedModelWithNoise.predictFcn(TI);
end
[HPR_max,idx] = max(HPR_static);
Qin_golden = Qtest(idx);

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
        % fprintf('\n Mechanistic HPR=%f, Inflow=%f, Gradient=%f',HPR_mechanistic(k),Qin_mechanistic(k-1),gradEst_mechanistic(k));
        
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


    %% Super Twisting - Extremun Seeking Control
    TI = table(GluIn(k), Qin_interpolated(k),'VariableNames',{'GluIn', 'Qin'});
    y2 = trainedModelWithNoise.predictFcn(TI);
    HPR = y2;
    
    if t(k)>=tic
        HPR_interpolated(k+1) = y2; 
        %fprintf('\n Interpolated HPR=%f, Inflow=%f, Gradient=%f',HPR_interpolated(k),Qin_interpolated(k-1),gradEst_interpolated(k));
        
        % Gradient approximation
        if t(k)>tic
            gradEst_interpolated(k+1) = (HPR_interpolated(k+1)-HPR_interpolated(k))/(Qin_interpolated(k)-Qin_interpolated(k-1));
        end
        
        % Optimization algorithm
        [~, QinOut2] = ode15s(@integratorST, [t(k) t(k+1)], Qin1_interpolated(k), optionsInt_interpolated,...
                            alphaST, gradEst_interpolated(k));
        Qin1_interpolated(k+1) = QinOut2(end);
        
        Qin_interpolated(k+1) = (lambdaST*sqrt(abs(gradEst_interpolated(k+1)))*sign(gradEst_interpolated(k+1))) + Qin1_interpolated(k+1);
    end


     %% PI
    TI_PI = table(GluIn(k), Qin_PI(k),'VariableNames',{'GluIn', 'Qin'});
    HPR_PI_temp = trainedModelWithNoise.predictFcn(TI_PI);

    if t(k) >= tic
        HPR_PI(k+1) = HPR_PI_temp;

        % Extremal PI error
        HPR_ref = HPR_max;
        e = HPR_ref - HPR_PI_temp;
        e_int = e_int + e * Dt;

        Qin_PI(k+1) =  Kp*e + Ki*e_int;

%         fprintf('\n PID-extremal HPR_max=%f, HPR=%f, Qin=%f, new Qin=%f, ...
%                 HPR_max, HPR_PID(k), Qin_PID(k), Qin_PID(k+1));
    end


    %% Classical ESC

    % Planta
    TI_ESC = table(GluIn(k), Qin_ESC(k), ...
                   'VariableNames', {'GluIn','Qin'});
    HPR_ESC_temp = trainedModelWithNoise.predictFcn(TI_ESC);

    if t(k)>=tic
        % Real output
        yC = HPR_ESC_temp;
        [~, eta] = ode15s(@HighPassFilterESC_ODE, ...
                          [t(k) t(k+1)], y1C, options, yC, wh);
        y1C = eta(end,1);
        HPR_ESC(k+1) = y1C;
        cita = (yC - y1C) * A * sin(w*t(k+1));

        [~, Qin_est] = ode15s(@IntegratorESC_ODE, ...
                              [t(k) t(k+1)], y2C, options, cita, kI);
        y2C = Qin_est(end,1);
        Qin_ESC(k+1) = y2C + A*sin(w*t(k+1));

%         % Saturation
%         Qin_ESC(k+1) = max(Qmin, min(Qmax, Qin_ESC(k+1)));
    end

%     fprintf('\n Classical ESC | HPR=%.4f | Qin=%.4f | New Qin=%.4f', ...
%             HPR_ESC(k), Qin_ESC(k), Qin_ESC(k+1));


    %% ST-ESC (proposed)
    TI_ST_ESC = table(GluIn(k), Qin_STESC(k),'VariableNames',{'GluIn', 'Qin'});
    HPR_SP_ESC_temp = trainedModelWithNoise.predictFcn(TI_ST_ESC);
    if t(k)>=tic
        HPR_STESC(k+1) = HPR_SP_ESC_temp;

        % Gradient approximation
        if t(k)>tic
            gradST(k+1) = (HPR_STESC(k+1)-HPR_STESC(k))/(Qin_STESC(k)-Qin_STESC(k-1));
        end
        
        % Optimization algorithm
        [~,qout] = ode15s(@integratorST,[t(k) t(k+1)],Qin1(k),options,alphaST,gradST(k));
        Qin1(k+1) = qout(end);
        Qin_STESC(k+1) = Qin1(k+1) + lambdaST*sqrt(abs(gradST(k)))*sign(gradST(k));
    end

end

% %% Plot results
% figure;
% plot(t, Qin_mechanistic, 'LineWidth', 2, 'MarkerSize', 6);
% hold on
% plot(t, Qin_interpolated, 'LineWidth', 2, 'MarkerSize', 6);
% xlabel('Time [d]', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Qin [L/d]', 'FontSize', 12, 'FontWeight', 'bold');
% %title('Input flow rate')
% legend('Qin to Mechanistic Model', ...
%        'Qin to Model without Cumulated Points', 'FontSize', 10)
% set(gca, 'FontSize', 12, 'LineWidth', 1.5)
% grid on
% box on
% 
% figure;
% plot(t, HPR_mechanistic, 'LineWidth', 2, 'MarkerSize', 6)
% hold on
% plot(t, HPR_interpolated, 'LineWidth', 2, 'MarkerSize', 6)
% xlabel('Time [d]', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('HPR [g[H_2]/Ld]', 'FontSize', 12, 'FontWeight', 'bold');
% %title('Hydrogen Production Rate')
% legend('HPR to Mechanistic Model', ...
%        'HPR to Model without Cumulated Points', 'FontSize', 10)
% set(gca, 'FontSize', 12, 'LineWidth', 1.5)
% grid on
% box on
% 
% figure;
% plot(t, gradEst_mechanistic, 'LineWidth', 2, 'MarkerSize', 6);
% hold on
% plot(t, gradEst_interpolated, 'LineWidth', 2, 'MarkerSize', 6);
% xlabel('Time [d]', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('\nabla HPR', 'FontSize', 12, 'FontWeight', 'bold');
% %title('Hydrogen Production Rate Gradient')
% legend('Gradient Mechanistic Model', ...
%        'Gradient Model without Cumulated Points', ...
%        'FontSize', 10)
% set(gca, 'FontSize', 12, 'LineWidth', 1.5)
% grid on
% box on


%% Plot HPR
figure; hold on; grid on
plot(t,HPR_PI,'c','LineWidth',2)
plot(t,HPR_ESC,'g','LineWidth',2)
plot(t,HPR_STESC,'b','LineWidth',2)
yline(HPR_max,'r--','LineWidth',2)
%yline(HPR_max,'r--','Golden optimum','LineWidth',2)
legend('PI','Classical ESC','ST-ESC','GSS')
%legend('PI','Classical ESC','ST-ESC','Golden Search')
xlabel('Time [d]')
ylabel('HPR [g[H_2]/Ld]')
set(gca,'FontSize',12)


%% Plot Qin
figure; hold on; grid on
plot(t,Qin_PI,'c','LineWidth',2)
plot(t,Qin_ESC,'g','LineWidth',2)
plot(t,Qin_STESC,'b','LineWidth',2)
yline(Qin_golden,'r--','LineWidth',2)
%yline(Qin_golden,'r--','Golden optimum','LineWidth',2)
legend('PI','Classical ESC','ST-ESC','GSS')
%legend('PI','Classical ESC','ST-ESC','Golden Search')
xlabel('Time [d]')
ylabel('Qin [L/d]')
set(gca,'FontSize',12)
