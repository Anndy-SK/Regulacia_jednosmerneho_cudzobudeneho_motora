%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% Zadanie:  Riadenie rychlosti jednosmerneho cudzobudeneho motora     %
% Autor:    Bc. Andrej Klein                                          %
% Datum:    16.4.2023                                                 %
% Skola:    Technicka Univerzita v Kosiciach, FEI, KEM                %
% Predmet:  Servopohony 2022/2023                                     %
% Prilohy:  Ku skriptu su vyhotovene tri simulacie v Simulinku        %
% Autor navrhu tohto skriptu, postupu a simulacii:  Bc. Andrej Klein  %
%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
clc, clear, close all, format compact; format short g

%% Uloha 1. Stavove veliciny, rychlost JCBM
% Parametre JCBM:
Ua = 420;
nN = 1410;
PN = 19500;
MN = 132;
Jm = 0.29;
Jc = Jm*8;
IN = 52;
Ra = 0.522;
La = 0.0081;
KTM = 100;
TTM = 5;

%% zakladne vzorce JCBM:
wN = (2*pi*nN)/60;
cfi = (Ua-Ra*IN)/wN;
Ui = cfi*wN;
Mz = MN;

%% stavove matice JCBM:
A = [0 (cfi/Jc);-(cfi/La) -(Ra/La)];
b = [0;1/La];
e = [-1/Jc;0];
cT = [1 0];

%% Uloha 2. Stavovy regulator s integracnym clenom
% SCHP
syms s r1 r2 Ki
I = eye(3);          % jednotkova matica 3.radu
rT = [r1 r2];
A_2 = [A-b*rT b; -Ki*cT 0];
SCHP = det(s*I - A_2);
SCHP = coeffs(SCHP,s);

% jednotlive koeficienty SCHP:
schp_f0 = SCHP(1,1)
schp_f1 = SCHP(1,2)
schp_f2 = SCHP(1,3)
schp_f3 = SCHP(1,4)

%% ZCHP
% zadanie: prekmit OS = 30%, doba regulacie tr = 1s
OS = 30; tr = 1;
% tlmenie d vyjadrene zo vzorca podla tabulky (po uprave):
d = sqrt((log(OS/100)^2)/(pi^2+log(OS/100)^2))
alpha = 5;  % nasobok

%% ZCHP pre 2% rozmedzie:
w0_2 = 4/(d*tr)
s1_2 = -d*w0_2 + (w0_2*sqrt(1-d^2))*i
s2_2 = -d*w0_2 - (w0_2*sqrt(1-d^2))*i
s3_2 = alpha*(-d*w0_2)
% zaokruhlime poly na 1. desatinne miesto
%s1_2 = round(s1_2,1)
%s2_2 = round(s2_2,1)
%s3_2 = round(s3_2,1)

ZCHP_2 = (s - s1_2)*(s - s2_2)*(s - s3_2);
ZCHP_2 = sym2poly(ZCHP_2)
zchp_f0_2 = ZCHP_2(1,4);
zchp_f1_2 = ZCHP_2(1,3);
zchp_f2_2 = ZCHP_2(1,2);
zchp_f3_2 = ZCHP_2(1,1);

r1_2 = solve(zchp_f1_2 == schp_f1);
r2_2 = solve(zchp_f2_2 == schp_f2);
Ki_2 = solve(zchp_f0_2 == schp_f0);

r1_2 = double(r1_2)
r2_2 = double(r2_2)
Ki_2 = double(Ki_2)

%% ZCHP pre 5% rozmedzie:
w0_5 = 1/(d*tr)*(3-0.5*log(1-d^2))
s1_5 = -d*w0_5 + (w0_5*sqrt(1-d^2))*i
s2_5 = -d*w0_5 - (w0_5*sqrt(1-d^2))*i
s3_5 = alpha*(-d*w0_5)
% zaokruhlime poly na 1. desatinne miesto
%s1_5 = round(s1_5,1)
%s2_5 = round(s2_5,1)
%s3_5 = round(s3_5,1)

ZCHP_5 = (s - s1_5)*(s - s2_5)*(s - s3_5);
ZCHP_5 = sym2poly(ZCHP_5)
zchp_f0_5 = ZCHP_5(1,4);
zchp_f1_5 = ZCHP_5(1,3);
zchp_f2_5 = ZCHP_5(1,2);
zchp_f3_5 = ZCHP_5(1,1);

r1_5 = solve(zchp_f1_5 == schp_f1);
r2_5 = solve(zchp_f2_5 == schp_f2);
Ki_5 = solve(zchp_f0_5 == schp_f0);

r1_5 = double(r1_5)
r2_5 = double(r2_5)
Ki_5 = double(Ki_5)

%% Uloha 3. Luenbergerov pozorovatel
% Riaditelnost a pozorovatelnost systemu:
QR = [b A*b];
if(det(QR) ~= 0)
    disp('System je riaditelny.');
else disp('System neni riaditelny.');
end
QP = [cT;cT*A]
if(det(QP) ~= 0)
    disp('System je pozorovatelny.');
else disp('System neni pozorovatelny.');
end
QP_inv = inv(QP)
qp = QP_inv(:,2)        % posledny stpec matice QP

%% VLASTNE HODNOTY SYSTEMU
disp('Vlastne hodnoty systemu:')
vlastne_hodnoty = eig(A)
vh1=vlastne_hodnoty(1)
vh2=vlastne_hodnoty(2)

%% KORENE POZOROVATELA:
disp('Korene pozorovatela:')
s1p = real(vh1)-5+i*imag(vh1)
s2p = real(vh2)-5+i*imag(vh2)

%% ZELANY POLYNOM POZOROVATELA
disp('Zelany polynom pozorovatela:')
zchr_p = poly([s1p s2p])
zchr_p_f2 = zchr_p(1)
zchr_p_f1 = zchr_p(2)
zchr_p_f0 = zchr_p(3)

%% URCENIE KOEFICIENTOV SPATNEJ VAZBY POZOROVATLA:
disp('Urcenie koeficientov spatnej vazby pozorovatela:')
I=eye(2)
h = [zchr_p_f0*I + zchr_p_f1*A + A^2]*qp
h1 = double(h(1))
h2 = double(h(2))