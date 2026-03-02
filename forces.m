% hello :)

%%%%%%%%% section for defining functions %%%%%%%%%%

function [u]=calcVelocity(H, w, k, h, z, t)
    u = H/2*w*((cosh(k*(h+z))/sinh(k*h)))*cos(w*t);
end

function [a]=calcAccel(H, w, k, h, z, t)
    a = H/2*w^2*((cosh(k*(h+z))/sinh(k*h)))*sin(w*t);
end

function [L]=setWidth(z)
    L = 10; % m
end

function [dF_drag]=calcDeltaDrag(H, w, k, h, z, t, C_D, dz)
rho = 1000; % kg/m^3
    u = calcVelocity(H, w, k, h, z+0.5*dz, t);
    z_1 = z + h;
    L = setWidth(z_1);
    A = dz.*L;
    dF_drag = 0.5 * rho .* u .* abs(u) .* C_D .* A;
end

function [dF_in]=calcDeltaInertia(H, w, k, h, z, t, C_M, dz)
    rho = 1000; % kg/m^3
    acc = calcAccel(H, w, k, h, z+0.5*dz, t);
    z_1 = z + h;
    V = dz .* setWidth(z_1) .* setWidth(z_1);
    dF_in = rho .* C_M .* V .* acc;
end

function [F_drag]=calcTotalDrag(H, w, k, h, t, C_D, dz, H_str)
    sum = 0;
    for z = 0:dz:H_str
        dF = calcDeltaDrag(H, w, k, h, z, t, C_D, dz);
        sum = sum + dF;
    end
    F_drag = sum;
end

function [F_iner]=calcTotalIner(H, w, k, h, t, C_M, dz, H_str)
    sum = 0;
    for z = 0:dz:H_str
        dF = calcDeltaInertia(H, w, k, h, z, t, C_M, dz);
        sum = sum + dF;
    end
    F_iner = sum;
end

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


%%%%%%%%% section for testing code %%%%%%%%%%

dt = 0.1;
t_end = 100;
H = 0.5; % meters
T = 10;
h = 1.5; % meters
z = -1; % meters

k = dispersionNR(pi, h, T, 0.01);
w = 2*pi / T

t_series = 0:dt:t_end;

u_series = calcVelocity(H, w, k, h, z, t_series);

% test the x-velocity
figure(1)
plot(t_series, u_series)
xlabel('Time (sec)')
ylabel('x-velocity (m/s)')
title('x-velocity at z = -1 m')

% test the drag force
C_D = 1.05;
dz = 0.01; % m
H_str = 1; 

drag_series = calcTotalDrag(H, w, k, h, t_series, C_D, dz, H_str);

figure(2)
plot(t_series, drag_series)
xlabel('Time (sec)')
ylabel('Drag (N)')
title('Drag force on structure')

% test the inertia force
C_M = 0.5;

iner_series = calcTotalIner(H, w, k, h, t_series, C_M, dz, H_str);

figure(3)
plot(t_series, iner_series)
xlabel('Time (sec)')
ylabel('Inertia (N)')
title('Inertia force on structure')

% test total force
figure(4)
total_force = drag_series+iner_series;
plot(t_series, total_force)
xlabel('Time (sec)')
ylabel('Total force (N)')
title('Total force on structure')