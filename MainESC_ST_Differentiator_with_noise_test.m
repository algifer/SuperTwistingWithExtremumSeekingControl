%% Simulation of the LIPATA dark fermentation reactor 
%% for hidrogen production (real process) with differentiator

clc; close all; clear;

SAVEFIGS    = 1;

%% Load the model parameters and experimental data

% Mecanistic Model
parametersLIPATA2014;
dataLIPATA2014;


%% Time definition
t0              = 0;   % d
tf              = 120;  %100;    % d
Dt              = 4;   % h

% Simulation time
t               = t0:Dt/(24):tf; %Sample time of 0.167 is each 4 hour
%t = t0:0.167:tf;
lt              = length(t);


%% Initial conditions
Glu0        = GluVec(1);            % g/L
Ace0        = AceVec(1);            % g/L
Pro0        = ProVec(1);            % g/L
Bu0         = BuVec(1);             % g/L
EtOH0       = EtOHVec(1);           % g/L
X0          = XVec(1);              % g/L

qCO2Vgas0   = qCO2VgasVec(1);       % L/d
qH2Vgas0    = qH2VgasVec(1);        % L/d

qCO2gas0    = ((Patm-pVapH2O)*qCO2Vgas0)/(R*TAmb*V);                   % mol/Ld
qH2gas0     = ((Patm-pVapH2O)*qH2Vgas0)/(R*TAmb*V);                    % mol/Ld

pCO2gas0    = (qCO2Vgas0/(qCO2Vgas0+qH2Vgas0)) * (Patm-pVapH2O);       % bar
CO2gas0     = pCO2gas0/(R*TReac);                                      % mol/L
CO20        = (qCO2gas0/KlaCO2) + (KHCO2*pCO2gas0);                    % mol/L

pH2gas0     = (qH2Vgas0/(qCO2Vgas0+qH2Vgas0)) * (Patm-pVapH2O);        % bar
H2gas0      = (pH2gas0/(R*TReac))*H2MM;                                % g/L
H20         = (qH2gas0/KlaH2) + (KHH2*pH2gas0*H2MM);                   % g/L

y0          = [Glu0; Ace0; Pro0; Bu0; EtOH0; X0; CO20; H20; CO2gas0; H2gas0];

%% ODEs solution
HPR         = zeros(1,lt);
Qin         = ones(1,lt);
Qin1        = ones(1,lt); % Input flowrate (L/d)
GluIn       = 25*ones(1, lt);
%GluIn_ode = linspace(5,30,lt); %20.0;         % g/L

y           = zeros(size(y0,1),lt);
y(:,1)      = y0;

%HPR(1)      = 0;

gradEst     = zeros(1,lt);
gradEst2    = zeros(1,lt);

% Optimization beginning
tic         = 5;                % d

% Extremum seeking parameters
alphaST     = 0.4; %0.0725;
lambdaST    = 0.091; %0.1;

flag        = 1;
e           = ones(4,1);       %Initial conditions for STA differentiator
options     = odeset('RelTol',1e-3,'AbsTol',1e-4);
optionsInt  = odeset('RelTol',1e-3,'AbsTol',1e-3);

[T, cOut]   = ode23(@H2ReactorLIPATA2RODE,t0:Dt/(24):tic, y(:,1),...
                      options, flag, GluIn(1), Qin(1), V, Vgas, K, H2MM, muMax, KGlu, KHCO2, KHH2, KlaCO2, KlaH2, Patm, pVapH2O, R, TReac, TAmb);
y(1:size(y0,1),1:length(T)) =     cOut.';
%HPR         = 100*y(10,:).*Qin./V;


GluIn = zeros(1, lt);
Qin_ref = zeros(1, lt);
HPR_ref = zeros(1, lt);
%% Noise and robustness metrics
HPR_meas    = zeros(1,lt);   % Measured HPR with noise
HPR_error   = zeros(1,lt);   % Error w.r.t. reference


% for k=1:lt
%     if((t(k)>=0) && (t(k)<5))
%         GluIn(k) = 0;
%         Qin_ref(k) = 0;
%         HPR_ref(k) = 0; 
%     elseif((t(k)>=5) && (t(k)<25))
%         GluIn(k) = 5;
%         Qin_ref(k) = 4.13;
%         HPR_ref(k) = 4.8; 
%     elseif((t(k)>=25) && (t(k)<50))
%         GluIn(k) = 20;
%         Qin_ref(k) = 4.74;
%         HPR_ref(k) = 24.6; 
%     elseif((t(k)>=50) && (t(k)<75))
%         GluIn(k) = 25;
%         Qin_ref(k) = 4.82;
%         HPR_ref(k) = 31.6; 
%     elseif((t(k)>=75) && (t(k)<95))
%         GluIn(k) = 10;
%         Qin_ref(k) = 4.6;
%         HPR_ref(k) = 11.4;
%     else
%         GluIn(k) = 15;
%         Qin_ref(k) = 4.66;
%         HPR_ref(k) = 17.8;
%     end
% end

for k=1:lt
    if t(k)>=0

        GluIn(k) = 25;
        Qin_ref(k) = 4.82;
        HPR_ref(k) = 31.6; 
    end
end

% --- Setup noise ---
noise_level = 0.05;                 % Noise level 1% (0.01)
noise_update_interval = 1;          % Each 1 second is updated the noise
rng(0);                            % Set the seed for reproducibility
noise_val = 1;                      % Initial value for noise

for k=1:lt-1
    if t(k)>=tic
       
       % Update the noise only when the base time changes
        if mod(t(k), noise_update_interval) == 0
            noise_val = 1 + noise_level * (2*rand - 1);  % new random value
        end
      
       [T, cOut]   = ode23(@H2ReactorLIPATA2RODEDiffReal,t(k):0.001:t(k+1), ..., 
                            [y(:,k); e], options, flag, GluIn(k), Qin(k), ...,
                            V, Vgas, K, H2MM, muMax, KGlu, KHCO2, KHH2, ..., 
                            KlaCO2, KlaH2, Patm, pVapH2O, R, TReac, TAmb, ...,
                            noise_val, al, be, ga, p, q, k1, k2);
        e           = cOut(end,end-3:end)';  % Save final conditions of the differentiator
        y(:,k+1)    = cOut(end,1:10)';       % Save final conditions of the state
        HPR(k)      = cOut(end,end-3);
        
        % Measurement noise (1% relative noise)
        HPR_meas(k) = HPR(k) * (1 + noise_level * (2*rand - 1));
        fprintf('\n Productivity=%f, Inflow=%f, Gradient=%f',HPR(k),Qin(k-1),gradEst(k));

        
        % Differentiator
        if t(k)>tic
            gradEst(k+1)   = cOut(end,end-1)./cOut(end,end);
        end

        % Optimization algorithm
            [~, QinOut]     = ode23tb(@integratorST, [t(k) t(k+1)], Qin1(k), ...
                optionsInt,alphaST, gradEst(k));
            Qin1(k+1)       = QinOut(end);
            Qin(k+1)        = (lambdaST*sqrt(abs(gradEst(k+1)))*sign(gradEst(k+1))) + Qin1(k);

    end
end

%% Quantitative noise robustness metrics

% Define steady-state window (last 20% of simulation)
ss_idx = round(0.8*lt):lt-1;

% Steady-state variance of HPR under noise
HPR_var = var(HPR_meas(ss_idx));

% RMS tracking error
HPR_rms_error = sqrt(mean((HPR_meas(ss_idx) - HPR_ref(ss_idx)).^2));

% Time to convergence (within 5% of reference)
tol = 0.05;
conv_idx = find(abs(HPR_meas - HPR_ref) ./ max(HPR_ref,1e-3) < tol, 1);

if isempty(conv_idx)
    t_conv = NaN;
else
    t_conv = t(conv_idx);
end

fprintf('\n--- Noise robustness metrics ---\n');
fprintf('Steady-state variance of HPR: %.4f\n', HPR_var);
fprintf('RMS tracking error: %.4f\n', HPR_rms_error);
fprintf('Time to convergence: %.2f days\n', t_conv);



%% Plot results
figure;
plot(t,GluIn,'LineWidth',2)
ylabel('GluIn [g/L]', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time [d]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter','tex');
set(gca, 'FontSize', 14, 'LineWidth', 1.5);
axis([0 max(t) 0 30])
grid on;
box on;

figure;
plot(t, Qin_ref, '--r', 'LineWidth',2)
hold on
plot(t,Qin,'LineWidth',2)
legend('Optimal','Q_{in}','Location','SouthEast')
ylabel('Qin [L/d]', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Time [d]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter','tex')
set(gca, 'FontSize', 14, 'LineWidth', 1.5);
grid on;
box on;

figure;
plot(t, HPR_ref, '--r', 'LineWidth',2)
hold on
plot(t(1:end-1), HPR(1:end-1),'LineWidth',2)
legend('Optimal','H_2 Production Rate','Location','SouthEast','FontSize',12)
ylabel('HPR [g[H_2]/Ld]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter','tex')
xlabel('Time [d]', 'FontSize', 14, 'FontWeight', 'bold')
set(gca, 'FontSize', 14, 'LineWidth', 1.5);
grid on;
box on;

figure;
plot(t,gradEst,'LineWidth',2)
ylabel('\nabla HPR', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Time [d]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter','tex')
set(gca, 'FontSize', 14, 'LineWidth', 1.5);
grid on;
box on;

figure;
% plot(t(1:end-1),HPR(1:end-1),'LineWidth',2)
% hold on
plot(t(1:end-1),HPR_meas(1:end-1),'--','LineWidth',2)
hold on
plot(t,HPR_ref,'k:','LineWidth',2)
%legend('True HPR','Measured HPR (1% noise)','Optimal','Location','SouthEast')
legend('Measured HPR (4% noise)','Optimal','Location','SouthEast')
ylabel('HPR [g H_2/Ld]')
xlabel('Time [d]')
grid on; box on;
