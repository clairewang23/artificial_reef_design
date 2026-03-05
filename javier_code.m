%% Morison horizontal force with vertical integration using generic A(z) and V'(z) vectors
% User provides vectors over height:
%   zVec  : elevations [m] (monotonic increasing), within [-h, 0]
%   bVec  : projected width normal to flow [m] at each z (array)
%   VpVec : displaced volume per unit height [m^2] at each z (Vp = dV/dz, array)
%
% Force eq. Morison:
%   dF_D = 0.5*rho*Cd*b(z)*dz*|u|u
%   dF_I = rho*Cm*Vp(z)*dz*u_dot
%
% Linear (Airy) wave, no current. z=0 at SWL, upward positive!!

clear; clc; close all;

%% -------------------- WAVE / WATER LEVEL --------------------
H   = 0.5;          % [m]
T   = 10;           % [s]
rho = 1025;         % [kg/m^3]
g   = 9.81;         % [m/s^2]
h   = 2.0;            % [m] water depth (seabed at z=-h)  <-- SET

nPeriods = 6;
dt       = 0.02;
% -------------------------------------------------------

%% -------------------- GEOMETRY AS VECTORS (GENERIC) --------------------
% Provide ANY shape by specifying zVec, bVec, VpVec (same length).
% Example below: arbitrary tapered element (edit/replace with your data).

reef_crest = -1.0;                        % [m] reef crest
zVec  = linspace(-h, -1.0, 5)';           % [m] elevations where geometry is defined (column)
bVec  = 1*1*ones(size(zVec));         % projected width normal to flow [m] at each z (array)
VpVec = 1^3*ones(size(zVec));           % displaced volume per unit height [m^2] at each z

% ----- Inputs for estability (DEFINE) -----
rho_s = 2400;   % material density [kg/m³] (concrete ~2400)
mu    = 0.60;   % base friction coefficient [-] (adjust according to substrate)
%Bbase = 1.0;    % Modify: width of the base in the direction of thrust [m] (resistance arm)
xCG   = 0.25;    % Modify: CG position relative to the leeard toe [m] 

rho_w = rho;    % usa la rho del agua del modelo

% Hydrodynamic Coefficients (to be modified)
Cd = 1.2;   
Cm = 2.0;   
% ------------------------------------------------------------------------


%% Wave number from dispersion
omega = 2*pi/T;
k = dispersionNR(0,h,T,0.001);

%% Time vector
tEnd = nPeriods*T;
t = 0:dt:tEnd;
phase = omega*t; % 1 x Nt

%% Kinematics u(z,t), u_dot(z,t) on zVec
Kz = cosh(k*(zVec + h)) / sinh(k*h);     % N x 1

u  = (H/2)*omega   * (Kz * cos(phase));  % N x Nt
ud = -(H/2)*omega^2 * (Kz * sin(phase)); % N x Nt

%% Vertical integration using trapezoidal rule (robust for non-uniform zVec)
% Drag integrand: 0.5*rho*Cd*b*|u|u
% Inertia integrand: rho*Cm*Vp*ud

% Precompute geometry weights (N x 1) to multiply each row
Gd = 0.5*rho .* (Cd .* bVec);   % N x 1
Gi = rho   .* (Cm .* VpVec);    % N x 1

% Build integrands (N x Nt)
fDrag = (Gd .* (abs(u).*u));       % N x Nt
fIner = (Gi .* ud);                % N x Nt

% Integrate over z for each t (trapz integrates along first dimension)
F_drag    = trapz(zVec, fDrag, 1); % 1 x Nt
F_inertia = trapz(zVec, fIner, 1); % 1 x Nt
F_total   = F_drag + F_inertia;

%% Plot: three forces in one figure
figure(1);
plot(t, F_drag,    'LineWidth', 1.2); hold on;
plot(t, F_inertia, 'LineWidth', 1.2);
plot(t, F_total,   'LineWidth', 1.6);
grid on; xlabel('t [s]'); ylabel('F [N]');
title('Morison: fuerzas horizontales integradas (geometría definida por vectores)');
legend('F_{drag}', 'F_{inercia}', 'F_{total}', 'Location', 'best');

%% Plot geometry vectors
figure(2); subplot(1,2,1)
plot(bVec, zVec, 'LineWidth', 1.2); grid on;
xlabel('b(z) [m]'); ylabel('z [m]'); title('Ancho proyectado b(z)');

figure(2); subplot(1,2,2)
plot(VpVec, zVec, 'LineWidth', 1.2); grid on;
xlabel('Vol (z)'); ylabel('z [m]'); title('Volumen por unidad de altura V''(z)');



%% =================== STABILITY: SLIDING + OVERTURNING ===================
% Requiere que ya existan:
% t (1xNt), zVec (Nx1), VpVec (Nx1), fDrag (NxNt), fIner (NxNt), F_total (1xNt)
% donde fDrag y fIner son integrandos [N/m] (antes de integrar en z), y F_total = ∫(fDrag+fIner) dz


% -----------------------------------------

%% Peso sumergido (W') usando V = ∫ V'(z) dz
Vtot = trapz(zVec, VpVec);                 % [m^3]
Wsub = (rho_s - rho_w) * 9.81 * Vtot;      % [N]  (W' del PDF)

%% SLIDING: FS_slide(t) = (mu*W') / |Fh(t)|
R_fric = mu * Wsub;                        % [N]
FS_slide_t = R_fric ./ max(abs(F_total), 1e-12);
FS_slide_min = min(FS_slide_t);

%% OVERTURNING: M_wave(t) = ∫ (z - z0) * dF(z,t)
% Momento respecto a la base (z0 = cota inferior del elemento)
z0 = zVec(1);                               % [m]
fTot = fDrag + fIner;                       % [N/m] integrando total
lever = (zVec - z0);                        % [m] brazo vertical
M_wave = trapz(zVec, (lever .* fTot), 1);   % [N·m], 1xNt

% Momento resistente: M_res = W' * e_res
% Brazo resistente genérico: desde el "toe" hasta la proyección del CG
% toe (borde) a +B/2 desde el centro. Si xCG>0 se acerca al toe.
e_res = xCG;                    % [m]

M_res = Wsub * e_res;                       % [N·m] (constante)

FS_overturn_t = M_res ./ max(abs(M_wave), 1e-12);
FS_overturn_min = min(FS_overturn_t);

%% Resumen
fprintf('\n=== STABILITY CHECKS ===\n');
fprintf('Vtot = %.4f m^3\n', Vtot);
fprintf('W''   = %.1f N\n', Wsub);
fprintf('Sliding:   min FS = %.3f  (FS(t)=mu*W''/|Fh(t)|)\n', FS_slide_min);
fprintf('Overturn:  min FS = %.3f  (FS(t)=M_res/|M_wave(t)|)\n', FS_overturn_min);

%% Plots
figure;
plot(t, FS_slide_t, 'LineWidth', 1.2); hold on;
plot(t, FS_overturn_t, 'LineWidth', 1.2);
yline(1.0, '--', 'LineWidth', 1.0);
grid on;
xlabel('t [s]'); ylabel('FS [-]');
title('Factores de seguridad: deslizamiento y vuelco');
legend('FS_{slide}', 'FS_{overturn}', 'FS=1', 'Location', 'best');

figure;
plot(t, F_total, 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('F_h [N]');
title('Fuerza horizontal total (Morison)');

figure;
plot(t, M_wave, 'LineWidth', 1.2); hold on;
yline(M_res, '--', 'LineWidth', 1.0);
grid on;
xlabel('t [s]'); ylabel('M [N·m]');
title('Momento de vuelco (onda) y momento resistente');
legend('M_{wave}(t)', 'M_{res}', 'Location', 'best');