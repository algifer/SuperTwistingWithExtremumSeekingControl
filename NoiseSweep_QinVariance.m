clc; clear;

noise_vec = 0:0.1:1.0;   % 0% a 100%
N = length(noise_vec);

Qin_var = zeros(N,1);
Qin_mean = zeros(N,1);

for i = 1:N
    nl = noise_vec(i);
    fprintf('Running noise level = %.0f %%\n', nl*100);

    [Qin, t] = Func_MainESC_ST_Differentiator_with_noise_test(nl);

    ss_idx = round(0.8*length(Qin)):length(Qin);
    Qin_var(i)  = var(Qin(ss_idx));
    Qin_mean(i) = mean(Qin(ss_idx));
end

Qin_var_rel = Qin_var / Qin_var(1);

%% ===== Results table =====
ResultsTable = table( ...
    noise_vec'*100, Qin_mean, Qin_var, Qin_var_rel, ...
    'VariableNames', {'Noise_percent','Mean_Qin','Var_Qin','Relative_Qin_variance'});

disp(ResultsTable)

