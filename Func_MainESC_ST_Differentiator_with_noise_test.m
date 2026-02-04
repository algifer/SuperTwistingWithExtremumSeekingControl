function [Qin_applied, t] = Func_MainESC_ST_Differentiator_with_noise_test(noise_level)
% ==========================================================
% ESC with ST differentiator – Noise injected in Qin channel
% ==========================================================

clc;

%% ================== SETUP ==================
parametersLIPATA2014;
dataLIPATA2014;

t0  = 0;
tf  = 120;
Dt  = 4/24;                % 4 h
t   = t0:Dt:tf;
lt  = length(t);

ticESC = 5;                % ESC activation time (days)

%% Initial conditions
Glu0 = GluVec(1);
Ace0 = AceVec(1);
Pro0 = ProVec(1);
Bu0  = BuVec(1);
EtOH0 = EtOHVec(1);
X0   = XVec(1);

qCO2Vgas0 = qCO2VgasVec(1);
qH2Vgas0  = qH2VgasVec(1);

qCO2gas0 = ((Patm-pVapH2O)*qCO2Vgas0)/(R*TAmb*V);
qH2gas0  = ((Patm-pVapH2O)*qH2Vgas0)/(R*TAmb*V);

pCO2gas0 = (qCO2Vgas0/(qCO2Vgas0+qH2Vgas0))*(Patm-pVapH2O);
CO2gas0  = pCO2gas0/(R*TReac);
CO20     = (qCO2gas0/KlaCO2)+(KHCO2*pCO2gas0);

pH2gas0 = (qH2Vgas0/(qCO2Vgas0+qH2Vgas0))*(Patm-pVapH2O);
H2gas0  = (pH2gas0/(R*TReac))*H2MM;
H20     = (qH2gas0/KlaH2)+(KHH2*pH2gas0*H2MM);

y0 = [Glu0;Ace0;Pro0;Bu0;EtOH0;X0;CO20;H20;CO2gas0;H2gas0];

%% ================== STORAGE ==================
y = zeros(10,lt);
y(:,1) = y0;

Qin        = zeros(1,lt);
Qin1       = zeros(1,lt);
Qin_applied= zeros(1,lt);
gradEst    = zeros(1,lt);
HPR        = zeros(1,lt);

GluIn = 25*ones(1,lt);

% Initial Qin
Qin(1)         = 4.8;
Qin1(1)        = 4.8;
Qin_applied(1)= 4.8;

%% ================== ESC PARAMETERS ==================
alphaST  = 0.4;
lambdaST = 0.091;

flag = 1;
e = ones(4,1);

options    = odeset('RelTol',1e-3,'AbsTol',1e-4);
optionsInt = odeset('RelTol',1e-3,'AbsTol',1e-3);

%% ================== NOISE SETUP ==================
rng(0);                            % reproducible
noise_update_interval = Dt;        % update once per sample
noise_val = 1;

%% ================== SIMULATION LOOP ==================
for k = 1:lt-1

    % --------- Noise update (affects Qin channel) ----------
    if t(k) >= ticESC
        if mod(t(k),noise_update_interval) < 1e-9
            noise_val = 1 + noise_level*(2*rand-1);
        end
    else
        noise_val = 1;
    end

    Qin_noisy = Qin(k)*noise_val;

    % --------- Reactor + differentiator ----------
    [~, cOut] = ode23( ...
        @H2ReactorLIPATA2RODEDiffReal2, ...
        t(k):0.001:t(k+1), ...
        [y(:,k); e], ...
        options, ...
        flag, GluIn(k), Qin_noisy, ...
        V, Vgas, K, H2MM, muMax, KGlu, KHCO2, KHH2, ...
        KlaCO2, KlaH2, Patm, pVapH2O, R, TReac, TAmb, ...
        noise_val, al, be, ga_ST, p, q, k1, k2);

    y(:,k+1) = cOut(end,1:10)';
    e        = cOut(end,end-3:end)';
    HPR(k)   = cOut(end,end-3);

    % --------- Gradient estimation ----------
    if t(k) > ticESC
        gradEst(k+1) = cOut(end,end-1) / cOut(end,end);
    else
        gradEst(k+1) = 0;
    end

    % --------- ESC update ----------
    if t(k) >= ticESC
        [~, QinOut] = ode23tb( ...
            @integratorST, ...
            [t(k) t(k+1)], Qin1(k), ...
            optionsInt, alphaST, gradEst(k));

        Qin1(k+1) = QinOut(end);
        Qin(k+1)  = Qin1(k+1) + ...
            lambdaST*sqrt(abs(gradEst(k+1))) * sign(gradEst(k+1));
    else
        Qin1(k+1) = Qin1(k);
        Qin(k+1)  = Qin(k);
    end

    Qin_applied(k+1) = Qin(k+1);
end

end

