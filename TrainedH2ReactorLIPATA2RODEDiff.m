%No esta terminado

function dy = TrainedH2ReactorLIPATA2RODEDiff(t, y, flag, u, J, al, be, ga, p, q, k1, k2)

%% Variables
% State variables
est     = y(1:2);                 % <----- Estimado de la senal a diferenciar
dest    = y(3:4);                 % <----- Estimado de las derivadas 

% Dilution rate
% u: Qin, J: HPR (Tasa de produccion de H_2)

errEst  = est - [J; u];
phi1    = flag*(al*power(norm(errEst,2),-p) + be + ga*power(norm(errEst,2),q))*errEst;
phi2    = flag*(al*(1-p)*power(norm(errEst,2),-p) + be + ga*(1+q)*power(norm(errEst,2),q))*phi1;


%% ODE system
dy  = [-k1*phi1 + dest;        %  <------ Estimado del estado ([J; u])
       -k2*phi2];              %  <------ Estimado de las derivadas

end