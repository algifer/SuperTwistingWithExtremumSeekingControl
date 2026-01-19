%% Set of parameters of the hydrogen production reactor LIPATA

%% Physicochemical constants
DifH2 = 4.65;           % m^2/s
DifCO2 = 1.98;          % m^2/s
KaAce = 1.51e-8;        % M
KaBu = 1.74e-5;         % M
kABi = 1e8;             % 1/d
KaCO2 = 4.94e-7;        % M
KaPro = 1.32e-5;        % M
KHCO2 = 0.0271;         % M/bar
KHH2 = 7.38e-4;         % M/bar
KlaCO2 = 2000;          % 1/d
KlaH2 = 3065;           % 1/d
Patm = 1.013;           % bar
pVapH2O = 0.0557;       % bar
R = 8.314e-2;           % bar/(K M)
TReac = 308;            % K
TAmb = 298;             % K
Vtot = 1.25;            % L
V = 0.9;                % L
Vgas = Vtot - V;        % L
HAM = 1.00794;          % g/M
H2MM = 2 * HAM;         % g/M
CO2MM = 44.01;          % g/M


% %% Steady state conditions at input
% GluIn = 10.0;         % g/L
% Qin = 4.32;            % L/d


%% Model parameters
% Pseudo-stoichiometric matrix

% K = [-1.0000   -1.0000
%          0    0.2118
%     0.0023    0.0038
%     0.3696    0.1874
%     0.0093    0.0159
%     0.0829    0.1068
%     0.0049    0.0087
%     0.0351         0];

% K = [-1.0000   -1.0000
%          0    0.3238
%     0.0205    0.0174
%     0.3028    0.2737
%     0.0242    0.0219
%     0.0934    0.1135
%     0.0072    0.0066
%     0.0381         0];

% K = [-1.0000   -1.0000
%          0    0.3033
%     0.0168    0.0142
%     0.3227    0.2403
%     0.0209    0.0285
%     0.1021    0.1018
%     0.0069    0.0072
%     0.0379         0];

K = [-1.0000   -1.0000
         0    0.3238
    0.0205    0.0174
    0.3028    0.2737
    0.0242    0.0219
    0.0829    0.1068
    0.0049    0.0087
    0.0351         0];

k11 = K(1,1);
k12 = K(1,2);
k21 = K(2,1);
k22 = K(2,2);
k31 = K(3,1);
k32 = K(3,2);
k41 = K(4,1);
k42 = K(4,2);
k51 = K(5,1);
k52 = K(5,2);
k61 = K(6,1);
k62 = K(6,2);
k71 = K(7,1);
k72 = K(7,2);
k81 = K(8,1);
k82 = K(8,2);

% Kinetic coefficients

% muMax1 = 35;                % g[Glu]/(g[X] d)
% KGlu1 = 0.18;               % g[Glu]/L
% muMax2 = 41;                % g[Glu]/(g[X] d)
% KGlu2 = 0.22;               % g[Glu]/L

muMax1 = 37.3197;           % g[Glu]/(g[X] d)
KGlu1 = 0.2896;             % g[Glu]/L
muMax2 = 27.2416;           % g[Glu]/(g[X] d)
KGlu2 = 0.2596;             % g[Glu]/L

% muMax1 = 25.9663;           % g[Glu]/(g[X] d)
% KGlu1 = 0.2578;             % g[Glu]/L

% muMax2 = 19.2578;           % g[Glu]/(g[X] d)
% KGlu2 = 0.2258;             % g[Glu]/L

muMax = [muMax1 muMax2];
KGlu = [KGlu1 KGlu2];

%Differeniator parameters
al      = 1; 
be      = 4;
ga      = 10;

p       = 0.35;
q       = 0.2;

k1      = 10;
k2      = 10;

% %Differeniator parameters
% al      = 1; 
% be      = 4;
% ga      = 5;
% 
% p       = 0.2;
% q       = 0.35;
% 
% k1      = 3;
% k2      = 3;
