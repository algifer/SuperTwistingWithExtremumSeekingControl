function dy = H2ReactorLIPATA2RODEDiffReal(t, y, flag, GluIn, Qin, V, Vgas, K, H2MM, muMax, KGlu, KHCO2, KHH2, KlaCO2, KlaH2, Patm, pVapH2O, R, TReac, TAmb, noise_val, al, be, ga_ST, p, q, k1, k2)

%% Global variables
global qGasCur

%% Variables
% State variables
Glu     = y(1,1);
Ace     = y(2,1);
Pro     = y(3,1);
Bu      = y(4,1);
EtOH    = y(5,1);
X       = y(6,1);
CO2     = y(7,1);
H2      = y(8,1);
CO2gas  = y(9,1);
H2gas   = y(10,1);

est     = y(11:12);                 % <----- Estimado de la senal a diferenciar
dest    = y(13:14);                 % <----- Estimado de las derivadas 

%% Gas phase
% Mass transfer
pGasCO2 = CO2gas*R*TReac;
pGasH2  = (H2gas*R*TReac)/H2MM;

rhoTCO2 = KlaCO2 * (CO2-(KHCO2*pGasCO2));
rhoTH2  = KlaH2 * (H2-(H2MM*KHH2*pGasH2));

qGasCur = ((R*TAmb)/(Patm-pVapH2O))*V*((rhoTH2/H2MM)+rhoTCO2);

% ODEs
dCO2gas = -((CO2gas*qGasCur)/Vgas) + (rhoTCO2*(V/Vgas));
dH2gas  = -((H2gas*qGasCur)/Vgas) + (rhoTH2*(V/Vgas));


% Input - Output (objective function)
D       = (Qin/V);  % u: Qin
J       = ((R*TAmb)/(Patm-pVapH2O))*((rhoTH2/H2MM));  % J: HPR

J = J * noise_val;

% % Add noise of 1%
% noise_level = 0.01; % 1%
% J = J + noise_level * J .* (2*rand - 1);
% J = J + noise_level * J .* randn(size(J)); 


errEst  = est - [J; Qin]; 
phi1    = flag*(al*power(norm(errEst,2),-p) + be + ga_ST*power(norm(errEst,2),q))*errEst;
phi2    = flag*(al*(1-p)*power(norm(errEst,2),-p) + be + ga_ST*(1+q)*power(norm(errEst,2),q))*phi1;

%% Liquid phase
% Reaction rate
r       = [(muMax(1)*Glu)/(KGlu(1)+Glu); (muMax(2)*Glu)/(KGlu(2)+Glu)] * X;

% ODEs
dGlu    = -(D*(Glu-GluIn)) + (K(1,:)*r);
dAce    = -(D*Ace) + (K(2,:)*r);
dPro    = -(D*Pro) + (K(3,:)*r);
dBu     = -(D*Bu) + (K(4,:)*r);
dEtOH   = -(D*EtOH) + (K(5,:)*r);
dX      = -(D*X) + (K(6,:)*r);
dCO2    = -(D*CO2) + (K(7,:)*r) - rhoTCO2;
dH2     = -(D*H2) + (K(8,:)*r) - rhoTH2;


%% ODE system
dy  = [[dGlu; dAce; dPro; dBu; dEtOH; dX; dCO2; dH2; dCO2gas; dH2gas];
       -k1*phi1 + dest;        %  <------ Inicia Estimador
       -k2*phi2];

end