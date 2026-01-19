function dy = H2ReactorLIPATA2RODE(t, y, flag, GluIn, Qin, V, Vgas, K, H2MM, muMax, KGlu, KHCO2, KHH2, KlaCO2, KlaH2, Patm, pVapH2O, R, TReac, TAmb)

%% Global variables
global qGasCur;

%% Variables
% State variables
Glu = y(1,1);
Ace = y(2,1);
Pro = y(3,1);
Bu = y(4,1);
EtOH = y(5,1);
X = y(6,1);
CO2 = y(7,1);
H2 = y(8,1);
CO2gas = y(9,1);
H2gas = y(10,1);

% Dilution rate
D = (Qin/V);

%% Gas phase
% Mass transfer
pGasCO2 = CO2gas*R*TReac;
pGasH2 = (H2gas*R*TReac)/H2MM;

rhoTCO2 = KlaCO2 * (CO2-(KHCO2*pGasCO2));
rhoTH2 = KlaH2 * (H2-(H2MM*KHH2*pGasH2));

qGasCur = ((R*TAmb)/(Patm-pVapH2O))*V*((rhoTH2/H2MM)+rhoTCO2);

% ODEs
dCO2gas = -((CO2gas*qGasCur)/Vgas) + (rhoTCO2*(V/Vgas));
dH2gas = -((H2gas*qGasCur)/Vgas) + (rhoTH2*(V/Vgas));


%% Liquid phase
% Reaction rate
r = [(muMax(1)*Glu)/(KGlu(1)+Glu); (muMax(2)*Glu)/(KGlu(2)+Glu)] * X;


% ODEs
dGlu = -(D*(Glu-GluIn)) + (K(1,:)*r);
dAce = -(D*Ace) + (K(2,:)*r);
dPro = -(D*Pro) + (K(3,:)*r);
dBu = -(D*Bu) + (K(4,:)*r);
dEtOH = -(D*EtOH) + (K(5,:)*r);
dX = -(D*X) + (K(6,:)*r);
dCO2 = -(D*CO2) + (K(7,:)*r) - rhoTCO2;
dH2 = -(D*H2) + (K(8,:)*r) - rhoTH2;


%% ODE system
dy = [dGlu; dAce; dPro; dBu; dEtOH; dX; dCO2; dH2; dCO2gas; dH2gas];

end