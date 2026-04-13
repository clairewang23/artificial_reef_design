clear all
close all

% Testing the design of our structure design. 

%% Define elements of the structure and other given values

% Constants
g = 9.81;                   % m/s^2, gravity acceleration
rho_w = 1025;               % kg/m^3, density of water
rho_s = 2400;               % kg/m^3, density of the structure

% Wave properties
h = 1.5;                         % m, water depth
H = 0.5;                         % m, wave height
T = 10;                          % s, wave period

% Geometry of the design


% INSERT VECTORS HERE
n = 10;                         % number of slices
height_str = 0.5;
dz = height_str/n;
V_vec = [0.00854136, 0.00980724, 0.0078201,  0.00395748, 0.00229142, ...
    0.00424953, 0.00289229, 0.0018537,  0.00145812, 0.00156665];
A_vec = [0.022, 0.02498838, 0.02225588, 0.01574212, 0.01269764, ...
    0.01643214, 0.0142229,  0.0107709,  0.00941921, 0.00992575];

%% Dispersion relation

% solve dispersion relation for k
function [Ksol]=dispersionNR(k0,h,T,eror)
    %       [Ksol]=dispersionNR(k0,h,T,eror);
    % 
    % k0     : Valor inicial de la iteracion (m-1)
    % h      : Valor de la profundidad (m)
    % T      : Valor del periodo (s)
    % eror   : Valor del error admitido
    
    w=2*pi/T;
    %k=2*pi./L;
    g=9.81;
    
    f=(w^2-g*k0.*tanh(k0*h));
    fp=-g*tanh(k0*h)-h*g*k0*(1-tanh(k0*h)^2);
    x1=1e10;
    x=0;
    i=0;
    
    while abs(x1-x)>eror  
       x=x1;
       i=i+1;
       f=(w^2-g*x.*tanh(x*h));
       fp=-g*tanh(x*h)-h*g*x*(1-tanh(x*h)^2);
       x1=x-f/fp;
       %pause
    end
Ksol=x1;
end

% Solve dispersion equation for wave number
k = dispersionNR(pi, h, T, 0.01);
w = 2*pi / T;

%% Functions to solve for forces

% calculate velocity
function [u]=calcVelocity(H, w, h, k, z, t)
    u = H/2*w*((cosh(k*(h+z))/sinh(k*h)))*cos(w*t);
end

% calculate acceleration
function [a]=calcAccel(H, w, h, k, z, t)
    a = H/2*w^2*((cosh(k*(h+z))/sinh(k*h)))*sin(w*t);
end

% calculate drag force differential given a specific height
function [dF_drag]=calcDeltaDrag(H, w, k, h, z, t, C_D, dz, rho_w, A)
    u = calcVelocity(H, w, k, h, z+0.5*dz, t);
    dF_drag = 0.5 * rho_w .* u .* abs(u) .* C_D .* A;
end

% calculate inertia force differential given specific height
function [dF_in]=calcDeltaInertia(H, w, k, h, z, t, C_M, dz, rho_w, V)
    acc = calcAccel(H, w, k, h, z+0.5*dz, t);
    dF_in = rho_w .* C_M .* V .* acc;
end

% calculate total drag force
function [F_drag]=calcTotalDrag(H, w, k, h, t, C_D, dz, rho_w, A_vec, n)
    sum = 0;
    for i = 1:n
        z = -h + dz*i;
        dF = calcDeltaDrag(H, w, k, h, z, t, C_D, dz, rho_w, A_vec(i));
        sum = sum + dF;
    end
    F_drag = sum;
end

% calculate total inertia force
function [F_iner]=calcTotalIner(H, w, k, h, t, C_M, dz, rho_w, V_vec, n)
    sum = 0;
    for i = 1:n
        z = -h + dz*i;
        dF = calcDeltaInertia(H, w, k, h, z, t, C_M, dz, rho_w, V_vec(i));
        sum = sum + dF;
    end
    F_iner = sum;
end

% calculate moment caused by wave forces
function [M]=calcMoment(H, w, k, h, t, C_D, C_M, dz, rho_w, ...
            A_vec, V_vec, n)
    sum = 0;
    for i = 1:n
        z = -h + dz*i;
        deltaI = calcDeltaInertia(H, w, k, h, z, t, C_M, dz, rho_w, V_vec(i));
        deltaD = calcDeltaDrag(H, w, k, h, z, t, C_D, dz, rho_w, A_vec(i));
        moment_arm = dz*i + 0.5*dz;
        dF = (deltaI + deltaD) * moment_arm;
        sum = sum + dF;
    end
    M = sum;    
end

%% Calculations
dt = 0.1;
t_end = 25;

t_series = 0:dt:t_end;
z = -0.7; % meters
u_series = calcVelocity(H, w, h, k, z, t_series);

% uncertain coefficients
C_D = 1;
C_M = 2;

% test velocity
figure(1)
plot(t_series, u_series), hold on
xlabel('Time (sec)')
ylabel('x-velocity (m/s)')
title('x-velocity at top of structure (z=-0.7m)')

% calculate the drag force
drag_series = calcTotalDrag(H, w, k, h, t_series, C_D, dz, rho_w, A_vec, n);

figure(2)
plot(t_series, drag_series, 'b'), hold on

% calculate the inertia force
iner_series = calcTotalIner(H, w, k, h, t_series, C_M, dz, rho_w, V_vec, n);

plot(t_series, iner_series, 'r'), hold on

% calculate and show total force
total_force = drag_series+iner_series;
plot(t_series, total_force, 'g')
xlabel('Time (sec)')
ylabel('Force (N)')
title('Forces on structure')

% add legend
legend('Drag force', 'Inertia force', 'Total force')

% print maximum horizontal force experienced by structure
max_force = max(total_force);
fprintf('Maximum horizontal force: %f N', max_force)

% test moment
moment_series = calcMoment(H, w, k, h, t_series, C_D, C_M, dz, rho_w, ...
    A_vec, V_vec, n);
figure(3)
plot(t_series, moment_series)
xlabel('Time (sec)')
ylabel('Moment (N*m)')
title('Moment on structure from wave action')

% print maximum moment on structure from waves
max_moment = max(moment_series);
fprintf('Maximum moment: %f N*m', max_moment)