% Use Morison equation to compute forces on the reef

%% Defining functions

% calculate velocity
function [u]=calcVelocity(H, w, k, h, z, t)
    u = H/2*w*((cosh(k*(h+z))/sinh(k*h)))*cos(w*t);
end

% calculate acceleration
function [a]=calcAccel(H, w, k, h, z, t)
    a = H/2*w^2*((cosh(k*(h+z))/sinh(k*h)))*sin(w*t);
end

%% Define functions for different shapes
% width function for a cube
function [L]=setWidth(z)
    L = 1 - 0.3*z; % m
end

% width of a sliced cone 
function [W]=setConeWidth(z, R1, R2, h)
    % z = height above the bottom of the element
    % R1 = radius of the bottom of the cone
    % R2 = radius of the top of the cone
    % h = height of the cone

    % calculate the width of the cone
    W = 2 * R1 - ((R1 - R2)/h)*z;
end

function [A]=setConeArea(z, R1, R2, h)
    % z = height above the bottom of the element
    % R1 = radius of the bottom of the cone
    % R2 = radius of the top of the cone
    % h = height of the cone

    % calculate the area of the cone
    A = pi * R1^2 - ((R1 - R2)/h) * z * 2 * pi;
end

% calculate drag force differential given a specific height
function [dF_drag]=calcDeltaDrag(H, w, k, h, z, t, C_D, dz)
rho = 1000; % kg/m^3
    u = calcVelocity(H, w, k, h, z+0.5*dz, t);
    z_1 = z + h;
    L = setWidth(z_1);
    A = dz.*L;
    dF_drag = 0.5 * rho .* u .* abs(u) .* C_D .* A;
end

% calculate inertia force differential given specific height
function [dF_in]=calcDeltaInertia(H, w, k, h, z, t, C_M, dz)
    rho = 1000; % kg/m^3
    acc = calcAccel(H, w, k, h, z+0.5*dz, t);
    z_1 = z + h;
    V = dz .* setWidth(z_1) .* setWidth(z_1);
    dF_in = rho .* C_M .* V .* acc;
end

% calculate total drag force
function [F_drag]=calcTotalDrag(H, w, k, h, t, C_D, dz, H_str)
    sum = 0;
    for z = 0:dz:H_str
        dF = calcDeltaDrag(H, w, k, h, z, t, C_D, dz);
        sum = sum + dF;
    end
    F_drag = sum;
end

% calculate total inertia force
function [F_iner]=calcTotalIner(H, w, k, h, t, C_M, dz, H_str)
    sum = 0;
    for z = 0:dz:H_str
        dF = calcDeltaInertia(H, w, k, h, z, t, C_M, dz);
        sum = sum + dF;
    end
    F_iner = sum;
end

% calculate moment caused by wave forces
function [M]=calcMoment(H, w, k, h, t, C_D, C_M, dz, H_str)
    sum = 0;
    for z = 0:dz:H_str
        deltaI = calcDeltaInertia(H, w, k, h, z, t, C_M, dz);
        deltaD = calcDeltaDrag(H, w, k, h, z, t, C_D, dz);
        dF = (deltaI + deltaD) * z;
        sum = sum + dF;
    end
    M = sum;    
end

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


%% Testing code

dt = 0.1;
t_end = 100;
H = 0.5; % meters
T = 10;
h = 2; % meters
z = -2; % meters

k = dispersionNR(pi, h, T, 0.01);
w = 2*pi / T;

t_series = 0:dt:t_end;

u_series = calcVelocity(H, w, k, h, z, t_series);

% test the x-velocity
%plot(t_series, u_series), hold on
%xlabel('Time (sec)')
%ylabel('x-velocity (m/s)')
%title('x-velocity at z = -1 m')

% test the drag force
C_D = 1.2;
dz = 0.01; % m
H_str = 1; 

drag_series = calcTotalDrag(H, w, k, h, t_series, C_D, dz, H_str);

figure(1)
plot(t_series, drag_series, 'b'), hold on

% test the inertia force
C_M = 2;

iner_series = calcTotalIner(H, w, k, h, t_series, C_M, dz, H_str);

plot(t_series, iner_series, 'r'), hold on

% test total force
total_force = drag_series+iner_series;
plot(t_series, total_force, 'g')
xlabel('Time (sec)')
ylabel('Force (N)')
title('Forces on structure')

% add legend
legend('Drag force', 'Inertia force', 'Total force')

% test moment
moment_series = calcMoment(H, w, k, h, t_series, C_D, C_M, dz, H_str);
figure(2)
plot(t_series, moment_series)
xlabel('Time (sec)')
ylabel('Moment (N*m)')
title('Momemt on structure from wave action')

% print maximum horizontal force experienced by structure
max_force = max(total_force);
fprintf('Maximum horizontal force: %f N', max_force)