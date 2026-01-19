parametersLIPATA2014;

%% Data from the hydrogen production reactor LIPATA 2014

%       Days(d)	Gluin(g/L)	Qin(ml/min)	OLR(g/Ld)	Glu(g/L) EtOH(g/L)	Ace(g/L)	Pro(g/L)	Bu(g/L) Biomass(g/L)	CO2(%)	H2(%)   qCO2(ml/h)	qH2(ml/h)	CODbalance(%)

data = [0.5     12.78	1.90	41.44	2.29	0.218	1.461	0.083	2.633	0.936	36.32	63.68	222.28	389.72	97.74
        ];
													
tVec = data(:,1)';                          % d
%Inputs
GluInVec = data(:,2)';                      % g/L
QinVec = data(:,3)' * (60*24/1000);         % L/d

GluVec = data(:,5)';                        % g/L
EtOHVec = data(:,6)';                       % g/L
AceVec = data(:,7)';                        % g/L
ProVec = data(:,8)';                        % g/L
BuVec = data(:,9)';                         % g/L
XVec = data(:,10)';                         % g/L

qCO2VgasVec = data(:,13)' * (24/1000);      % L/d
qH2VgasVec = data(:,14)' * (24/1000);       % L/d

HprVec = qH2VgasVec./V;
qGasVec = qCO2VgasVec + qH2VgasVec;         % L/d
perCO2gasVec = (qCO2VgasVec./qGasVec)*100;  % Percentage
perH2gasVec = (qH2VgasVec./qGasVec)*100;    % Percentage
CODBalance = data(:,15)';                   % Percentage

% Outputs
qCO2gas = ((Patm-pVapH2O)*qCO2VgasVec)/(R*TAmb*V)';                   % mol/Ld
qH2gas = ((Patm-pVapH2O)*qH2VgasVec)/(R*TAmb*V)';                     % mol/Ld

% QinSim = QinVec'
% GluInSim = GluInVec'
% qH2VgasSim = qH2VgasVec';
% HprSim = qH2VgasSim./V
% [QinInput, GluInput] = meshgrid(QinSim,GluInSim);
% xvec = QinSim;
% yvec = GluInSim;
% [Xa, Ya] = ndgrid(xvec, yvec);
% F = scatteredInterpolant(xvec, yvec, HprSim);
% Za = F(Xa, Ya);
% figure (100)
% surf(Xa, Ya, Za, 'edgecolor', 'none');
