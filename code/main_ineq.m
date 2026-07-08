clear;
clc;
close all

set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14);
set(0,'DefaultFigureWindowStyle', 'docked'); 
set(0,'defaulttextInterpreter','latex');
rng('default');

%% 0. Parametri del DC-DC Boost Converter 
% Parametri del sistema
T = 0.65e-3;    % [s]
L = 4.2e-3;     % [H]
C = 2200e-6;    % [F]
Res = 85;       % [Ohm]
Vin = 15;       % [V]

% Stati e ingressi principali
x_ref = [0.389; -16];   u_ref = 0.516;    % Punto di equilibrio
x_0 = [2.3; -13.5];                         % Punto iniziale

% Parametri MPC 
T_sim = 60;     % Numero di step simulati
N = 20;         % Orizzonte di predizione
Q = eye(2);     % Matrice costo Q
R = 1;          % Costo R

% Vincoli fisici assoluti
iL_min = -1.6;      iL_max = 2.4;       
vo_min = -19;       vo_max = -13;       
dc_min = 0.3;       dc_max = 0.75;  

%% 1. Invariant set per il sistema controllato

% Sistema linearizzato al discreto valutato in (x_ref, u_ref)
A = [1,                    (T/L)*(1 - u_ref);
    -(T/C)*(1 - u_ref),    1 - (T/(Res*C))];

B = [-(T/L)*(x_ref(2) - Vin);
     (T/C)*x_ref(1)];

% Vincoli su stato e ingresso (nelle coordinate del sistema linearizzato)
Hx = [eye(2); -eye(2)];
hx = [iL_max - x_ref(1); 
      vo_max - x_ref(2); 
     -iL_min + x_ref(1); 
     -vo_min + x_ref(2)];

Hu = [1; -1];
hu = [dc_max - u_ref; 
     -dc_min + u_ref];

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

% Riferimento del sistema linearizzato
x_ref_lin = [0; 0];
u_ref_lin = 0;

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

%% 5. Plot temporali (in millisecondi)
figure(4);
set(gcf, 'Color', 'w'); 

% Parametri grafici uniformi
lw = 1.2;              
legFontSize = 8;      
grayCol = [0.6 0.6 0.6]; 

% Calcolo del vettore tempo in millisecondi
t_ms = (0:T_sim) * T * 1000; 
t_max_ms = t_ms(end);

% Sottografico per x1 (Corrente iL)
ax(1) = subplot(3,1,1);
p1 = plot(t_ms, x_log(1,:), 'LineWidth', lw, 'Color', 'b');
hold on;
p2 = yline(x_ref(1), '--', 'Color', grayCol, 'LineWidth', lw);
p3 = yline(iL_min, '--r', 'LineWidth', lw);
yline(iL_max, '--r', 'LineWidth', lw, 'HandleVisibility', 'off');

ylabel('$x_1: i_L$ [A]');
title('\textbf{Evoluzione temporale degli stati e dell''ingresso}');
grid on; xlim([0 t_max_ms]);
ylim([iL_min - 0.15*abs(iL_min), iL_max + 0.15*abs(iL_max)]);
legend([p1, p2, p3], 'Stato x_1', 'Equilibrio', 'Vincoli', ...
       'FontSize', legFontSize, 'Location', 'northeast');

% Sottografico per x2 (Tensione vo) 
ax(2) = subplot(3,1,2);
p1 = plot(t_ms, x_log(2,:), 'LineWidth', lw, 'Color', 'b');
hold on;
p2 = yline(x_ref(2), '--', 'Color', grayCol, 'LineWidth', lw);
p3 = yline(vo_min, '--r', 'LineWidth', lw);
yline(vo_max, '--r', 'LineWidth', lw, 'HandleVisibility', 'off');

ylabel('$x_2: v_o$ [V]');
grid on; xlim([0 t_max_ms]);
ylim([vo_min - 0.1*(abs(vo_max-vo_min)), vo_max + 0.1*(abs(vo_max-vo_min))]);
legend([p1, p2, p3], 'Stato x_2', 'Equilibrio', 'Vincoli', ...
       'FontSize', legFontSize, 'Location', 'northeast');

% Sottografico per u (Duty Cycle dc) 
ax(3) = subplot(3,1,3);
p1 = stairs(t_ms(1:end-1), u_log, 'LineWidth', lw, 'Color', 'b');
hold on;
p2 = yline(u_ref, '--', 'Color', grayCol, 'LineWidth', lw); 
p3 = yline(dc_min, '--r', 'LineWidth', lw);
yline(dc_max, '--r', 'LineWidth', lw, 'HandleVisibility', 'off');

xlabel('Tempo [ms]');
ylabel('$u: dc$ [-]');
grid on; xlim([0 t_max_ms]);
ylim([dc_min - 0.05, dc_max + 0.05]);
legend([p1, p2, p3], 'Ingresso u', 'Equilibrio', 'Vincoli', ...
       'FontSize', legFontSize, 'Location', 'northeast');

