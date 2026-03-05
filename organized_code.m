% Testing the design of a particular structure. 

%% Define elements of the structure

% Constants
g = 9.81;                   % m/s^2, gravity acceleration
rho_w = 1025;               % kg/m^3, density of water
rho_s = 2400;               % kg/m^3, density of the structure

% Properties of the structure (cone)
R = 1/2;                          % m, radius of bottom of structure
r = 0.4/2;                        % m, radius of top of structure
h = 0.85;                         % m, height of the structure
V = (pi*h*(R^2 + R*r + r^2))/3;   % volume, m^3
W = (rho_s - rho_w)*g*V;          % weight of the structure

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





% Wave properties
