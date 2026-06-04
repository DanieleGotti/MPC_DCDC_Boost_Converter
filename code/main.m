clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
set(0,'defaulttextInterpreter','latex')
rng('default');

%%  0. Parametri del DC-DC Boost Converter (in SI units)
T = 0.65e-3;    % [s]
L = 4.2e-3;     % [H]
C = 2200e-6;    % [F]
R = 85;         % [Ohm]
Vin = 15;       % [V]

%%  1. Invariant set per il sistema controllato

% Riferimento FISICO del sistema (Equilibrio)
x_ref = [0.389; -16];
u_ref = 0.516e-3; 

% Sistema linearizzato al discreto valutato in (x_ref, u_ref)
A = [1,                         (T/L)*(1 - u_ref);
    -(T/C)*(1 - u_ref),          1 - (T/(R*C))];

B = [-(T/L)*(x_ref(2) - Vin);
     (T/C)*x_ref(1)];

% Vincoli fisici assoluti
iL_min = -1.6;      iL_max = 2.4;       
vo_min = -19;       vo_max = -13;       
dc_min = 0.3e-3;    dc_max = 0.75e-3;   

% Vincoli su stato e ingresso (NELLE COORDINATE DEL SISTEMA LINEARIZZATO!)
% Trasliamo i vincoli sottrando l'equilibrio
Hx = [eye(2); -eye(2)];
hx = [iL_max - x_ref(1); 
      vo_max - x_ref(2); 
     -iL_min + x_ref(1); 
     -vo_min + x_ref(2)];

Hu = [1; -1];
hu = [dc_max - u_ref; 
     -dc_min + u_ref];

% Matrici del costo quadratico
Q = eye(2);
R = 1;

% Computazione control invariant set