clear;
clc;
close all

set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14);
set(0,'DefaultFigureWindowStyle', 'docked'); 
set(0,'defaulttextInterpreter','latex');
rng('default');

%% 0. Parametri del DC-DC Boost Converter 
T = 0.65e-3;    % [s]
L = 4.2e-3;     % [H]
C = 2200e-6;    % [F]
Res = 85;       % [Ohm]
Vin = 15;       % [V]

%% 1. Invariant set per il sistema controllato

% Riferimento fisico del sistema (equilibrio)
x_ref = [0.389; -16];
u_ref = 0.516; 

% Sistema linearizzato al discreto valutato in (x_ref, u_ref)
A = [1,                    (T/L)*(1 - u_ref);
    -(T/C)*(1 - u_ref),    1 - (T/(Res*C))];

B = [-(T/L)*(x_ref(2) - Vin);
     (T/C)*x_ref(1)];

% Vincoli fisici assoluti
iL_min = -1.6;      iL_max = 2.4;       
vo_min = -19;       vo_max = -13;       
dc_min = 0.3;       dc_max = 0.75;  

% Vincoli su stato e ingresso (nelle coordinate del sistema linearizzato)
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
[CIS_H, CIS_h] = cis(A, B, [0; 0], 0, Hx, hx, Hu, hu, Q, R);
CIS = Polyhedron(CIS_H, CIS_h);

figure(1)
CIS.plot('Alpha', 0.6);
title('\textbf{Control invariant set sistema linearizzato}');
xlabel('$\delta i_L$ [A]');
ylabel('$\delta v_o$ [V]');
grid on;

%% 2. N-step-controllable set dell'invariant set

% Orizzonte di predizione
N = 10; 

% Calcolo dell'insieme controllabile a N passi
[Np_steps_H, Np_steps_h] = controllable_set(Hx, hx, Hu, hu, CIS_H, CIS_h, A, B, N);
Np_step_set = Polyhedron('A', Np_steps_H, 'b', Np_steps_h);

figure(2)
Np_step_set.plot('Alpha', 0);
title('\textbf{CIS e N-step-controllable set}');
xlabel('$\delta i_L$ [A]');
ylabel('$\delta v_o$ [V]');
hold on;
CIS.plot('Alpha', 0.6);
grid on;

%% 3. Design MPC e simulazione

% Numero di step simulati
T_sim = 60;

% Riferimento del sistema linearizzato
x_ref_lin = [0; 0];
u_ref_lin = 0;

% Stato iniziale
x_0 = [2.3; -13.5];

mpc = mpc_ingredients_ineq(A, B, Hx, hx, Hu, hu, CIS_H, CIS_h, x_ref_lin, u_ref_lin, Q, R, N);

% Log stati e ingresso sistema
x_log = zeros(2,T_sim+1);
u_log = zeros(1,T_sim);
flags = zeros(1,T_sim);

x_log(:,1) = x_0;

for tt = 1:T_sim

    % Stato sistema linearizzato
    x_lin = x_log(:,tt) - x_ref;
    x_lin_shifted = x_lin - x_ref_lin;

    % Impostazioni MPC relative alla condizione iniziale
    f = mpc.f_base * x_lin_shifted;
    b_ineq = mpc.b_ineq_base - mpc.b_ineq_x0_factor * x_lin_shifted;
    
    % Risoluzione problema di ottimizzazione
    [delta_u_seq, ~, exitflag] = quadprog(mpc.F, f, mpc.A_ineq, b_ineq);

    flags(tt) = exitflag;

    if exitflag ~= 1
        warning('Quadprog fallito allo step %d. (Exitflag: %d).', tt, exitflag);
    end

    % Azione di controllo fisica
    u_log(tt) = u_ref + delta_u_seq(1);
    
    % Risposta del sistema
    x_log(:,tt+1) = boost_converter(x_log(:,tt), u_log(tt), T, L, C, Res, Vin);

end

%% 4. Plot dei risultati

% Traslazione del CIS e dell'N-step set nelle coordinate originali
CIS_shifted = CIS + x_ref;
Np_step_set_shifted = Np_step_set + x_ref;

figure(3);
Np_step_set_shifted.plot('Alpha',0);
title('\textbf{Traiettoria del sistema}');
xlabel('$i_L$ [A]');
ylabel('$v_o$ [V]');
hold on;
CIS_shifted.plot('Alpha', 0.6);
hold on;
plot(x_log(1,:),x_log(2,:),'Color',[0 0 0.5]);
scatter(x_log(1,:),x_log(2,:),'cyan');
grid on;