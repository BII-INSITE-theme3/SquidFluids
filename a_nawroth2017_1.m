%% Reconstruction of the two cylinder flow field

% Author: Stephen Williams

% History: Code made (07/06/23)

% 

%% Preamble

 clear all
% close all

%% Constansts

% System parameters
Xsi = 200; NgridX = 100;
Ysi = 200; NgridY = 100;
a = 90; % Cyclinder diameter (um)
Nstok = 100; % Number of stokeslets per cylinder
eps = 3; % Regularisation parameter (um)
cent1x = -a; cent1y = 0; % Coordinates of cyclinder 1 center
cent2x = a; cent2y = 0; % Coordinates of cyclinder 2 center

% Fluid parameters
mu = 1.08e-3; % Fluid viscocity (Pa.s)
Oo4PiMu = 1/(4*pi*mu);

%% Preallocation

f0 = (8*pi)/(1-2*log(a))* [0,1]; % (8*pi)/(1-2*log(a))* [1,0]; % Change this later
theta = linspace(0,2*pi,Nstok); % Angles at which stokeslets lie on cyclinder (rad)
stokesCoords1 = zeros(2,Nstok);
stokesCoords2 = zeros(2,Nstok);
stokesCoords1 = [cent1x + (a*0.95/2).*sin(theta);cent1y + (a*0.95/2).*cos(theta)]; % Left cylinder coords
stokesCoords2 = [cent2x + (a*0.95/2).*sin(theta);cent2y + (a*0.95/2).*cos(theta)]; % Right cylinder coords

X = linspace(-Xsi,Xsi,NgridX); % X grid points to calculate u
Y = linspace(-Ysi,Ysi,NgridY); % Y grid points to calculate u
Ux = zeros(NgridX,NgridY);
Uy = zeros(NgridX,NgridY);

%% Main loop

for nx = 1:NgridX % Loop on x coords
    for ny = 1:NgridY % Loop on y coords

        pos = [X(nx),Y(ny)]; % Get current position
        
        for s = 1:Nstok

            Lstok = stokesCoords1(:,s); % coords for left stokeslet
            Rstok = stokesCoords2(:,s); % coords for right stokeslet

            rL = norm(Lstok - pos);
            rR = norm(Rstok - pos);

            % Get contributions from left cylinder
            Ux(nx,ny) = Ux(nx,ny) + ...
                     -f0(1)*Oo4PiMu*( log(sqrt(rL^2 + eps^2) + eps ) - eps*(sqrt(rL^2 + eps^2) - 2*eps) )/((sqrt(rL^2 + eps^2) - 2*eps)*sqrt(rL^2 + eps^2)) + ...
                     Oo4PiMu * (f0(1) * (pos(1) - Lstok(1)) + f0(2) * (pos(2) - Lstok(2)))* f0(1) * (pos(1) - Lstok(1)) * (sqrt(rL^2 + eps^2) - 2*eps)/((sqrt(rL^2 + eps^2) - 2*eps)^2*sqrt(rL^2 + eps^2));
            Uy(nx,ny) = Uy(nx,ny) + ...
                     -f0(2)*Oo4PiMu*( log(sqrt(rL^2 + eps^2) + eps ) - eps*(sqrt(rL^2 + eps^2) - 2*eps) )/((sqrt(rL^2 + eps^2) - 2*eps)*sqrt(rL^2 + eps^2)) + ...
                     Oo4PiMu * (f0(1) * (pos(1) - Lstok(1)) + f0(2) * (pos(2) - Lstok(2)))* f0(2) * (pos(2) - Lstok(1)) * (sqrt(rL^2 + eps^2) - 2*eps)/((sqrt(rL^2 + eps^2) - 2*eps)^2*sqrt(rL^2 + eps^2));
             % Get contributions from right cylinder
            Ux(nx,ny) = Ux(nx,ny) + ...
                     -f0(1)*Oo4PiMu*( log(sqrt(rR^2 + eps^2) + eps ) - eps*(sqrt(rR^2 + eps^2) - 2*eps) )/((sqrt(rR^2 + eps^2) - 2*eps)*sqrt(rR^2 + eps^2)) + ...
                     Oo4PiMu * (f0(1) * (pos(1) - Rstok(1)) + f0(2) * (pos(2) - Rstok(2)))* f0(1) * (pos(1) - Rstok(1)) * (sqrt(rR^2 + eps^2) - 2*eps)/((sqrt(rR^2 + eps^2) - 2*eps)^2*sqrt(rR^2 + eps^2));
            Uy(nx,ny) = Uy(nx,ny) + ...
                     -f0(2)*Oo4PiMu*( log(sqrt(rR^2 + eps^2) + eps ) - eps*(sqrt(rR^2 + eps^2) - 2*eps) )/((sqrt(rR^2 + eps^2) - 2*eps)*sqrt(rR^2 + eps^2)) + ...
                     Oo4PiMu * (f0(1) * (pos(1) - Rstok(1)) + f0(2) * (pos(2) - Rstok(2)))* f0(2) * (pos(2) - Rstok(1)) * (sqrt(rR^2 + eps^2) - 2*eps)/((sqrt(rR^2 + eps^2) - 2*eps)^2*sqrt(rR^2 + eps^2));

        end

    end
end

%% Plotting tools

%% Functions 
