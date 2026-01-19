function evaluate_and_save_points(points, bounds)
    % Define simulation time
    t0 = 0;   % d
    tf = 30;  % d
    Dt = 4;   % d
    t = t0:Dt/(24):tf;
    lt = length(t);
    
    % Get inputs
    GluIn_DOE = points(:, 1);
    Qin_DOE = points(:, 2);
    
    % Number of tests
    testNumber = 3;
    points_count = size(points, 1);

    % Initialize variables with noise
    GluIn_with_noise = zeros(testNumber, points_count);
    Qin_with_noise = zeros(testNumber, points_count);
    t1 = 1:points_count;
    
    % Generate noise
    for jj = 1:testNumber
        GluIn_noise = GluIn_DOE .* (1 + 3 * 2 * (randn(points_count, 1) - 0.5) / 100);
        Qin_noise = Qin_DOE .* (1 + 1 * 2 * (randn(points_count, 1) - 0.5) / 100);
        for kk = 1:points_count
            GluIn_with_noise(jj, kk) = GluIn_noise(kk);
            Qin_with_noise(jj, kk) = Qin_noise(kk);
        end
    end
 
    % Load model parameters and data
    parametersLIPATA2014;
    dataLIPATA2014;
    
    % Data size
    [rData, cData] = size(data);
    datIni = 1;
    datFin = rData;
    
    % Global variables
    global qGasCur;
    
    % Condiciones iniciales
    Glu0 = GluVec(1);    % g/L
    Ace0 = AceVec(1);    % g/L
    Pro0 = ProVec(1);    % g/L
    Bu0 = BuVec(1);      % g/L
    EtOH0 = EtOHVec(1);  % g/L
    X0 = XVec(1);        % g/L
    
    qCO2Vgas0 = qCO2VgasVec(1); % L/d
    qH2Vgas0 = qH2VgasVec(1);   % L/d
    
    qCO2gas0 = ((Patm - pVapH2O) * qCO2Vgas0) / (R * TAmb * V); % mol/Ld
    qH2gas0 = ((Patm - pVapH2O) * qH2Vgas0) / (R * TAmb * V);   % mol/Ld
    
    pCO2gas0 = (qCO2Vgas0 / (qCO2Vgas0 + qH2Vgas0)) * (Patm - pVapH2O); % bar
    CO2gas0 = pCO2gas0 / (R * TReac);                                  % mol/L
    CO20 = (qCO2gas0 / KlaCO2) + (KHCO2 * pCO2gas0);                    % mol/L
    
    pH2gas0 = (qH2Vgas0 / (qCO2Vgas0 + qH2Vgas0)) * (Patm - pVapH2O);  % bar
    H2gas0 = (pH2gas0 / (R * TReac)) * H2MM;                           % g/L
    H20 = (qH2gas0 / KlaH2) + (KHH2 * pH2gas0 * H2MM);                 % g/L
    
    y0 = [Glu0; Ace0; Pro0; Bu0; EtOH0; X0; CO20; H20; CO2gas0; H2gas0];
    
    % ODEs solution
    Hpr_with_noise = zeros(testNumber, points_count);
    
    % Integrator Options
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    
    for j = 1:testNumber
        GluIn = GluIn_with_noise(j, :);
        for k = 1:points_count
            % Call IVP solver
            [~, yOut] = ode15s('H2ReactorLIPATA2RODE', t, y0, options, GluIn(k), Qin_DOE(k), V, Vgas, K, ...
                H2MM, muMax, KGlu, KHCO2, KHH2, KlaCO2, KlaH2, Patm, pVapH2O, R, TReac, TAmb);
            
            % Simulation results
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
            Hpr_value = QH2gas/V;
            Hpr_noise = Hpr_value*(1 + 3*2*(randn(1)-0.5)/100);
            Hpr_with_noise(j,k) = Hpr_noise;
        end
    end
    
    % Averages of HPr and GluIn
    GluIn_mean_vect = mean(GluIn_with_noise, 1);
    Hpr_mean_vect = mean(Hpr_with_noise, 1);

    % Create table with data
    T = table(GluIn_mean_vect', Qin_DOE, Hpr_mean_vect', 'VariableNames', {'GluIn', 'Qin', 'Hpr'});
    
    % Save table to an Excel file
    fileName = sprintf('IO_regression_By_DOE_with_Noise_%d_cumulated.xlsx', points_count);
    writetable(T, fileName, 'Sheet', 'Hoja1', 'Range', 'A1');
    
    % Plot the results
    figure;
    plot3(GluIn_mean_vect, Qin_DOE, Hpr_mean_vect, 'o')
    ylabel('Input flow rate [L/d]');
    xlabel('Input glucose [g/L]');
    zlabel('HPr [L(H_2)L^{-1}d^{-1}]');
end