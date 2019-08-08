%% Detailing settings for simulation
rho=1.225;
c = 343;
% Simulating Scale Model from jensens paper. Details found in chapter 13.

% vacuum permeability
u0=4*pi*1e-7;
% relative recoil permeability
ur =  1.12; %  1.05-1.2.
% Br = Effective magnetic remanence i.e. the intersection between the recoil line and the ordinate (B-field) axis
Br = 0.41; % 0.4-0.42 T
% Fm=Br*Dm/u0; -> Br=Fm*u0/Dm;
% A= Cross-sectional area of permanent magnet and air gap when they are assumed to be the same 
Am=144/10^6; % 12mm * 12mm
Ag=72/10^6; % 12mm * 6mm
A=Ag;
% D = Length of one air gap when the armature is in its resting position [
D =  1.9/1000;
% lm is the magnet length. 
lm = 10/1000;
Dm= Ag/Am*lm/ur;
% Effective D
Deff=D+Dm;
% Number of coil windings
N = 1755;
% Electromotive Force
Fm =Br*Dm/u0;
% Electrical DC resistance of armature, R_L.
R_e= 10.7;


% length of armature, chapter 13 joe jensen
l_a=90/1000;
% height (or thickness) of armature, chapter 13 joe jensen
h_a=1.35/1000;
% width of armature, chapter 13 joe jensen
w_a=6/1000;
% volume of moving part of armature
v1=l_a*h_a*w_a;
% volume of the added extrusion, c.a. 2cm long and looks like a cube
v2=h_a^2*20/1000;
% rho_a is the armature mass density, assuming magnetic steel???
rho_a=7700;
% Equivalent armature mass Meq
Meq=33/140*rho_a*(v1+v2);
% M = mechanical moving mass representing the armature and any other relevant contributions
% M = Meq;
% fn= first fundamental frequency.
% fn=1/(2*pi)*sqrt(ka(0)*1000/Meq);
M=k_a(0)/(128*2*pi)^2;

% r=0.11;

% settings defines what approximations and nonlinearities are used in the
% 1=joe, 2 = klippel, 3=adjustedklippel
settings = 1;
settingsStringHelper={'Joe','Klippel','adjustedKlippel'};
settingsString=strcat('SimSettings(',settingsStringHelper(settings),')');

% simulation setting for which nonlinearityto include
% T_em and T_me are closely related and combined to T. same with
% T_emd,T_med to T_d
% nonLinSettings = [k_o, T, T_d, L];
% Every combination of nonLinSettings:
nonLinSettings=dec2bin(0:15)-'0';
stringComboHelper = {'k_o-','T-','Td-','L-'};

% [forward euler  or midpoint, sumthin]
simSettings=[true,false];
