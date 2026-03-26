clear;
close all;
clc;

%% def variabili e valori costanti
% In questa sezione si preparano:
% - variabili simboliche (per derivare Jacobiane e linearizzare in modo analitico)
% - costanti numeriche del problema (parametri fisici dati in tabella)
% x1 = theta (angolo), x2 = omega (velocità angolare)
% u  = Cm (coppia motore, ingresso)
% J_i, psi_i = coefficienti/fasi delle armoniche del momento d'inerzia J(theta)
% beta = attrito viscoso, k = costante elastica, theta_e = posizione di equilibrio
%definizione variabili simboliche 

syms x1 x2 u J0 J1 J2 J3 J4 psi1 psi2 psi3 psi4 beta k theta_e real

%valori numerici tratti dalle specifiche
k_val=50;
beta_val=1;
J0_val=1.7;
J1_val=0.1;
J2_val=0.7;
J3_val=0.09;
J4_val=0.02;
psi1_val=-0.04;
psi2_val=2.9;
psi3_val=2.8;
psi4_val=-2.6;

% valore dell'angolo di equilibrio assegnato
theta_e_val=pi/3;

% Questi due estremi vengono usati più avanti (nel punto 5 della consegna)
% per esplorare un range di condizioni iniziali del modello non lineare.
% Qui vengono solo definiti.
theta_e_val_max = theta_e_val + 28;
theta_e_val_min = theta_e_val - 29; 

%% PUNTO 1: modello
% Punto 1 della consegna:
% - riscrivere il modello in forma di stato xdot=f(x,u), y=h(x,u)
% - trovare (xe,ue) di equilibrio
% - linearizzare: matrici A,B,C,D (Jacobiane valutate in equilibrio)

%definizione J(theta)
J_x1 = J0 + J1*cos(x1 + psi1) + J2*cos(2*x1 + psi2) + J3*cos(3*x1 + psi3) + J4*cos(4*x1 + psi4);

%funzioni di stato e di uscita
% Scelta dello stato coerente con la consegna/relazione:
% x1 := theta, x2 := omega
% Equazioni:
%  x1_dot = omega = x2
%  x2_dot = omega_dot = (u - beta*omega - k*theta)/J(theta)
f1 = x2;
f2 = (u - beta*x2 - k*x1) / J_x1;

% vettore della funzione di stato f(x,u)
f = [f1; f2];
% uscita misurata: y = theta = x1
h = x1;

%%calcolo matrici simboliche
% Linearizzazione: δxdot = A δx + B δu ; δy = C δx + D δu
% dove:
% A = ∂f/∂x |(xe,ue)
% B = ∂f/∂u |(xe,ue)
% C = ∂h/∂x |(xe,ue)
% D = ∂h/∂u |(xe,ue)

%calcolo matrice jacobiana A
% Derivate parziali di f1 e f2 rispetto a x1 e x2
df1_dx1 = diff(f1, x1);
df1_dx2 = diff(f1, x2);
df2_dx1 = diff(f2, x1);
df2_dx2 = diff(f2, x2);

% Assemblaggio della jacobiana A (2x2)
A_sym = [df1_dx1, df1_dx2;
         df2_dx1, df2_dx2];

%calcolo matrice jacobiana B
% Derivate di f rispetto all'ingresso u
df1_du = diff(f1, u);
df2_du = diff(f2, u);

% Assemblaggio della jacobiana B: è un vettore colonna (2x1)
B_sym = [df1_du; df2_du];

%calcolo matrice jacobiana C
% Derivate dell'uscita h rispetto allo stato: qui h=x1 quindi C=[1 0]
dh_dx1 = diff(h, x1);
dh_dx2 = diff(h, x2);

% Assemblaggio della jacobiana C: è un vettore riga (1x2)
C_sym = [dh_dx1, dh_dx2];

% Calcolo della matrice jacobiana D 
% Derivata dell'uscita rispetto all'ingresso: qui h non dipende da u quindi D=0
dh_du = diff(h, u);

% Assemblaggio della jacobiana D: è uno scalare
D_sym = dh_du;

%definizione punto di equilibrio
% Equilibrio: f(xe,ue)=0.
% Condizioni:
%  x1e = theta_e (dato)
%  x2e = 0 (velocità a regime nulla)
%  Da x2_dot=0 -> (ue - beta*0 - k*theta_e)/J(theta_e)=0 -> ue = k*theta_e

x1_e = theta_e;
x2_e = 0;
x_e = [x1_e; x2_e];
u_e = k * theta_e;

%definizione uscita di equilibrio
% Valori numerici dell'equilibrio (servono per stampare e creare sistema numerico)
y_e_val = theta_e_val;
u_e_val = k_val * theta_e_val;
x_e_val = [theta_e_val; 0];

% Questi due vettori vengono usati più avanti (nel punto 5 della consegna)
% per esplorare un range di condizioni iniziali del modello non lineare.
% Qui vengono solo definiti.
x_e_val_max = [theta_e_val_max; 0]; 
x_e_val_min = [theta_e_val_min; 0]; 

%calcolo delle matrici nel punto di equilibrio
% Sostituisco (x1,x2,u) con (x1e,x2e,ue) nelle jacobiane simboliche
A_eq = subs(A_sym, [x1, x2, u], [x1_e, x2_e, u_e]);
B_eq = subs(B_sym, [x1, x2, u], [x1_e, x2_e, u_e]);
C_eq = subs(C_sym, [x1, x2, u], [x1_e, x2_e, u_e]);
D_eq = subs(D_sym, [x1, x2, u], [x1_e, x2_e, u_e]);

%calcoli numerici
% Calcolo J(theta_e)
J_e_val = J0_val + J1_val*cos(theta_e_val + psi1_val) +   J2_val*cos(2*theta_e_val + psi2_val) +   J3_val*cos(3*theta_e_val + psi3_val) + J4_val*cos(4*theta_e_val + psi4_val);
fprintf('J_e = %.6f\n\n', J_e_val);


% Qui sostituisco i parametri fisici (J0..J4, psi1..psi4, beta,k,theta_e)
% ottenendo A e B numeriche.
A_num = double(subs(A_eq, {J0, J1, J2, J3, J4, psi1, psi2, psi3, psi4, beta, k, theta_e}, {J0_val, J1_val, J2_val, J3_val, J4_val, psi1_val, psi2_val, psi3_val, psi4_val, beta_val, k_val, theta_e_val}));
B_num = double(subs(B_eq, {J0, J1, J2, J3, J4, psi1, psi2, psi3, psi4, beta, k, theta_e}, {J0_val, J1_val, J2_val, J3_val, J4_val, psi1_val, psi2_val, psi3_val, psi4_val, beta_val, k_val, theta_e_val}));

% C e D non dipendono dai parametri: si convertono direttamente in double

C_num = double(C_eq);
D_num = double(D_eq);

% Stampa delle matrici linearizzate numeriche
disp('Matrice A numerica:');
disp(A_num);

disp('Matrice B numerica:');
disp(B_num);

disp('Matrice C numerica:');
disp(C_num);

disp('Matrice D numerica:');
disp(D_num);

%% PUNTO 2: Funzione di trasferimento
% Punto 2 della consegna:
% ricavare G(s) tra ingresso δu e uscita δy a partire dal modello linearizzato.

%solo per visualizzione, range di pulsazioni per i diagrammi di Bode
omega_plot_min = 1e-4;
omega_plot_max = 1e7;


%sistema linearizzato
% Creo il sistema in spazio di stato continuo (ss) con matrici numeriche
sys_lin = ss(A_num, B_num, C_num, D_num);

%funzione di trasferimento G
G = tf(sys_lin);
G %per stampare sul terminale il valore della fdt G (per controllo)

%calcolo dei poli e degli zeri
% Poli: autovalori del sistema e poli della fdt. Zeri: zeri della fdt.
p = pole(G);
disp('Poli della G');
disp(p);

z = zero(G);
disp('Zeri della G');
disp(z);

% Mappa poli-zeri: controllo visivo
figure(1);
grid on;
pzmap(G);

%Diagramma di bode
% Valuta modulo e fase nel range di pulsazioni definito sopra
figure(2);
bode(G, {omega_plot_min,omega_plot_max});

%% PUNTO 3: 
% Punto 3 della consegna:
% Mappatura specifiche -> vincoli in frequenza (zone proibite sul diagramma di Bode di L)
% Qui si impostano le specifiche e si traducono in regioni vietate (patch) su modulo/fase.

% 1) errore a regime
% limite massimo ammesso sull'errore a regime a gradino
e_reg_spec = 0.04;

% 2) margine di fase
% requisito minimo di robustezza (in gradi)
M_f_spec = 30;

% 3) tempo di assestamento
% Dal progetto: S% <= 10% -> S_100_spec = 0.10 (cioè 10% in forma decimale)
S_100_spec = 0.10;

% Dalla formula che lega sovraelongazione e smorzamento xi:
%   S% = exp(-π ξ / sqrt(1-ξ^2))
% che si inverte ottenendo la "xi_star" usata poi per stimare un margine di fase equivalente.
xi_star = abs(log(S_100_spec)) / sqrt(pi^2 + (log(S_100_spec))^2);
disp("xi_star");
disp(xi_star);
% Qui si usa la formula: Mf_star = 100*xi_star
% e si trova il valore di M_f_spec più stringente possibile da specifiche

M_f_spec = max(xi_star*100, M_f_spec);

% 4) tempo di assestamento
T_a5_spec = 0.008;

% 5) specifica disturbo d(t) (bassa frequenza)
% attenuazione richiesta: almeno 30 dB
Ad = 30;  %dB
omega_d_min = omega_plot_min; %lower bound
omega_d_MAX = 0.5;

% 6) specifica rumore n(t) (alta frequenza)
% attenuazione richiesta: almeno 65 dB
An = 65;  %dB
omega_n_min = 5*1e4;
omega_n_MAX = 5*1e6;

figure(3);
hold on

%specifiche su d
% Zona proibita:
patch_d_x = [omega_d_min; omega_d_min; omega_d_MAX; omega_d_MAX];
patch_d_y = [-300 ; Ad; Ad; -300];
patch(patch_d_x, patch_d_y, 'r', 'FaceAlpha', 0.1,'EdgeAlpha',0, 'DisplayName', 'Specica disturbo d(t)');

%specifiche su n
% Zona proibita:
patch_n_x = [omega_n_min; omega_n_min; omega_n_MAX; omega_n_MAX];
patch_n_y = [-An ; 300; 300; -An];
patch(patch_n_x, patch_n_y, 'b', 'FaceAlpha', 0.1,'EdgeAlpha',0, 'DisplayName', 'Specica rumore n(t)');

%specifiche tempo di assestamento
% Vincolo su velocità: per Ta piccolo serve una ωc sufficientemente grande.
omega_Ta_min = omega_plot_min;    %lower bound
omega_Ta_MAX = 300/(M_f_spec*T_a5_spec);

% Zona proibita:
patch_Ta_x = [omega_Ta_min; omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX];
patch_Ta_y = [-300 ; 0; 0; -300];
patch(patch_Ta_x, patch_Ta_y, 'g', 'FaceAlpha', 0.1,'EdgeAlpha',0, 'DisplayName', 'Vincolo tempo di assestamento');

%Plot Bode con margini di stabilità
margin(G,{omega_plot_min,omega_plot_max});
grid on; zoom on;

%specifiche sovraelongazione (margine di fase)
% Vincolo sul Mf viene rappresentato come "zona proibita" in fase tra ωc_min e ωn_min.
% Intervallo scelto:
% - minimo: omega_c_min = omega_Ta_MAX (da vincolo di velocità)
% - massimo: omega_c_MAX = omega_n_min (prima della banda alta dove devi attenuare rumore)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = omega_n_min;

% Per avere Mf >= M_f_spec:
% la fase in attraversamento deve essere >= -180 + Mf_spec

phi_up = M_f_spec-180;
phi_low = -270;
% Zona proibita:
% si crea patch tra phi_low e phi_up
patch_Mf_x = [omega_c_min; omega_c_min; omega_c_MAX; omega_c_MAX];
patch_Mf_y = [phi_low; phi_up; phi_up; phi_low];
patch(patch_Mf_x, patch_Mf_y, 'y', 'FaceAlpha', 0.3,'EdgeAlpha',0, 'DisplayName', 'Vincolo margine di fase');

legend show


%% PUNTO 4: Regolatore Statico
% Punto 4 (inizio progetto del regolatore):
% Qui si sintetizza un guadagno statico Rs per soddisfare:
% - errore a regime (L(0) sufficientemente grande)
% - reiezione disturbo d(t) alla banda bassa (|L(jωd,max)| sufficientemente grande)

%valori presi dalle specifiche
% ampiezze massime di riferimento e disturbo (da consegna): W<=1, D<=1
W = 1;
D = 1;

% valore minimo prescritto per L(0)
mu_s_error = (D+W)/e_reg_spec;

% vincolo disturbo: attenuazione Ad -> richiede |L(jωd,max)| >= 10^(Ad/20)
mu_s_dist  = 10^(Ad/20);

% guadagno minimo del regolatore ottenuto come L(0)/G(0)
% G_0 = |G(j0)|, G_omega_d_MAX = |G(jωd,max)|
G_0 = abs(evalfr(G,0));
G_omega_d_MAX = abs(evalfr(G,j*omega_d_MAX));

% Rs deve rendere:
%   Rs*|G(0)| >= µs_error
%   Rs*|G(jωd,max)| >= µs_dist
% quindi Rs >= max(µs_error/|G(0)|, µs_dist/|G(jωd,max)|)
R_s = max(mu_s_error/G_0,mu_s_dist/G_omega_d_MAX); 

%% Sistema esteso
G_e = R_s*G;

figure(4);
hold on

% Ripeto le stesse "zone proibite" e confronto su Bode di Ge
%specifiche su d
patch_d_x = [omega_d_min; omega_d_min; omega_d_MAX; omega_d_MAX];
patch_d_y = [-300 ; Ad; Ad; -300];
patch(patch_d_x, patch_d_y, 'r', 'FaceAlpha', 0.1,'EdgeAlpha',0, 'DisplayName', 'Specica disturbo d(t)');

%specifiche su n
patch_n_x = [omega_n_min; omega_n_min; omega_n_MAX; omega_n_MAX];
patch_n_y = [-An ; 300; 300; -An];
patch(patch_n_x, patch_n_y, 'b', 'FaceAlpha', 0.1,'EdgeAlpha',0, 'DisplayName', 'Specica rumore n(t)');

%specifiche tempo di assestaamento
patch_Ta_x = [omega_Ta_min; omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX];
patch_Ta_y = [-300 ; 0; 0; -300];
patch(patch_Ta_x, patch_Ta_y, 'g', 'FaceAlpha', 0.1,'EdgeAlpha',0, 'DisplayName', 'Vincolo tempo di assestamento');

%plot Bode con margini di stabilità
margin(G_e,{omega_plot_min,omega_plot_max});
grid on; zoom on;

%specifiche sovraelongazione (margine di fase)
patch_Mf_x = [omega_c_min; omega_c_min; omega_c_MAX; omega_c_MAX];
patch_Mf_y = [phi_low; phi_up; phi_up; phi_low];
patch(patch_Mf_x, patch_Mf_y, 'y', 'FaceAlpha', 0.3,'EdgeAlpha',0, 'DisplayName', 'Vincolo margine di fase');

legend show

%% regolatore dinamico SCENARIO B
% Dopo Rs, siamo nello "Scenario B":
% - specifiche a bassa frequenza ok
% - mancano margine di fase e/o ωc per soddisfare Ta
% Quindi si progetta una rete anticipatrice + filtro passa-basso ad alta freq.

% Scelta "con margine": Mf_star = Mf_spec + 10°, ωc_star = 1.1*ωc_min
% per robustezza e per non essere troppo vicini ai vincoli.
Mf_star = M_f_spec + 10;
omega_c_star = 1.1*omega_c_min;  

% Valuto modulo e fase di Ge alla pulsazione desiderata ωc_star
mag_omega_c_star = abs(evalfr(G_e,j*omega_c_star));
arg_omega_c_star = rad2deg(angle(evalfr(G_e,j*omega_c_star)));

% Obiettivo rete:
% - Portare |Ge(jωc_star) * Rd(jωc_star)| = 1  (0 dB)  -> M_star = 1/|Ge|
% - Aggiungere fase per arrivare a Mf_star in attraversamento.
M_star = 1/mag_omega_c_star;

% Fase necessaria:
% Mf = 180 + fase(L(jωc))   -> per avere Mf_star imposto:
% fase(Rd) = Mf_star - 180 - fase(Ge)
phi_star = Mf_star - 180 - arg_omega_c_star;

% Formule di inversione per rete anticipatrice:
% Rd(s) = (1 + τ s) / (1 + α τ s), con 0<α<1 per una rete di anticipo
% Qui si ricavano τ e α imponendo modulo e fase in ωc_star.
tau = (M_star - cos(phi_star*pi/180))/(omega_c_star*sin(phi_star*pi/180));
alpha_tau = (cos(phi_star*pi/180) - 1/M_star)/(omega_c_star*sin(phi_star*pi/180));
alpha = alpha_tau / tau;

% Controllo di consistenza: parametri fisicamente sensati (positivi)
if min(tau,alpha) < 0
    fprintf('Errore: parametri rete anticipatrice negativi');
    return;
end

%% Diagramma di Bode con  sistema in anello chiuso con regolatore
% Creo il regolatore dinamico completo:
% - rete anticipatrice
% - filtro passa-basso R_high_frequency per attenuare rumore in alta frequenza
s = tf('s');

% polo ad alta frequenza (1/(1+s/ω_hf)) con ω_hf=9e4 rad/s
R_high_frequency = 1/(1 + s/9e4);

% Rd(s) = rete anticipatrice (per aumentare il margine di fase) x
% filtro passa-basso (per limitare il guadagno alle alte frequenze)
R_d = (1 + tau*s)/(1 + alpha * tau*s)*R_high_frequency;

% Regolatore totale: R(s)=Rs*Rd(s)
R = R_s*R_d;

% Funzione d'anello L(s)=R(s)G(s)
L = R*G;


figure(5);
hold on;

% specifiche su ampiezza (stesse patch già definite)
patch(patch_d_x, patch_d_y, 'r', 'FaceAlpha', 0.1,'EdgeAlpha',0, 'DisplayName', 'Specica disturbo d(t)');
patch(patch_n_x, patch_n_y, 'b', 'FaceAlpha', 0.1,'EdgeAlpha',0, 'DisplayName', 'Specica rumore n(t)');
patch(patch_Ta_x, patch_Ta_y, 'g', 'FaceAlpha', 0.1,'EdgeAlpha',0, 'DisplayName', 'Vincolo tempo di assestamento');

% plot Bode con margini di stabilità per L(s)
margin(L,{omega_plot_min,omega_plot_max});
grid on; zoom on;

% specifiche su fase (vincoli Mf)
patch(patch_Mf_x, patch_Mf_y, 'y', 'FaceAlpha', 0.3,'EdgeAlpha',0, 'DisplayName', 'Vincolo margine di fase');

legend show;

%% Test sul sistema linearizzato
% Punto 4 della consegna (test sul linearizzato):
% Si simulano le risposte in anello chiuso a:
% - riferimento w(t)
% - disturbo d(t) (sommatoria armoniche a bassa freq)
% - rumore di misura n(t) (sommatoria armoniche ad alta freq)
% usando funzioni di sensitività.

% funzioni di sensitività
F=L/(1+L); % -> sensitività complementare 
S=1/(1+L); % -> sensitività

% Si costruiscono vettori tempo con passi diversi:
% - tt_n: passo molto piccolo per catturare rumore ad alta frequenza (wn ~ 5e4)
% - tt_d: simulazione lunga e passo grosso per disturbo a bassa frequenza (wd ~ 0.1)
% - tt: simulazione "principale" per risposta complessiva, con passo piccolo

t_simulazione_n=1e-3;
passo_n=1e-7;
tt_n=[0:passo_n:t_simulazione_n];

t_simulazione_d=300;
passo_d=0.01;
tt_d=[0:passo_d:t_simulazione_d];

t_simulazione=0.05;
passo=passo_n;
tt=[0:passo:t_simulazione];

% costruisco il segnale di riferimento
% Gradino unitario w(t)=1
WW=1;
ww=WW*ones(length(tt),1);

% costruisco i disturbi
% d(t) = sum_{k=1..4} D sin(wd*k*t)  (bassa frequenza)
% n(t) = sum_{k=1..4} N sin(wn*k*t)  (alta frequenza)
DD=1; dd=0; wd=0.1; dd_plot=0;
NN=1; nn=0; wn=5e4; nn_plot=0;

for k=[1:1:4]
    dd=dd + DD*sin(wd*k*tt); 
    nn=nn + NN*sin(wn*k*tt); 
    dd_plot=dd_plot + DD*sin(wd*k*tt_d);
    nn_plot=nn_plot+NN*sin(wn*k*tt_n);  
end 

% calcolo l'uscita per il segnale di riferimento
% y_w(t) = F(s) * w(t)
y_w=lsim(F,ww,tt);

% calcolo le uscite per il segnale di riferimento e per i disturbi
% Disturbo d(t) entra "in uscita": y_d = S * d
y_d=lsim(S, dd, tt);
y_d_plot=lsim(S, dd_plot, tt_d);

% Rumore n(t) entra sul ramo di misura e si sottrae: contribuzione in uscita -F*n
y_n=lsim(-F, nn, tt);
y_n_plot=lsim(-F,nn_plot,tt_n);

% eseguo il plot per la risposta al segnale di riferimento
figure(6);
grid on, zoom on, hold on;
plot(tt, ww, 'm'); % riferimento ideale
plot(tt, y_w, 'b');  % risposta all'ingresso w tramite F
legend('ww','y_w');

% eseguo il plot per la risposta al disturbo di attuazione
figure(7);
hold on, grid on, zoom on
plot(tt_d,dd_plot,'m'); % disturbo d(t)
plot(tt_d,y_d_plot,'b');  % uscita dovuta a d(t) filtrata da S
legend('dd','y_d');

% eseguo il plot per la risposta al disturbo di misura
figure(8);
hold on, grid on, zoom on;
plot(tt_n,nn_plot,'m');  % rumore n(t)
plot(tt_n,y_n_plot,'b');  % uscita dovuta a n(t) filtrata da -F
legend('nn','y_n');

% eseguo il plot sul sistema linearizzato considerando sia il segnale di
% riferimento che i rumori
% Principio di sovrapposizione degli effetti (poichè siamo in un sistema LTI)
y_tot = y_w + y_d + y_n;
figure(9);
grid on, zoom on, hold on;
plot(tt,ww,'m');
plot(tt, y_tot,'b');

% LV è il valore di regime previsto per l'uscita con ingresso WW.
% evalfr(WW*F,0) = WW*F(0) = guadagno in continua della catena chiusa.
LV = evalfr(WW*F,0); 

% vincolo su sovraelongazione/errore rispetto al valore di regime
patch([0, 0, t_simulazione, t_simulazione],  [LV*(1+S_100_spec), LV*2, LV*2, LV*(1+S_100_spec)], 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.5);
% vincolo su tempo di assestamento a +/- 5% dopo T_a5_spec
patch([T_a5_spec,T_a5_spec, t_simulazione, t_simulazione], [LV*(1-0.05),0, 0, LV*(1-0.05)], 'g', 'FaceAlpha', 0.1, 'EdgeAlpha', 0.5);
patch([T_a5_spec,T_a5_spec, t_simulazione, t_simulazione],[LV*(1+0.05),LV*2, LV*2, LV*(1+0.05)],'g', 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1);
legend('ww','y_tot');