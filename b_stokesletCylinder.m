%% Script to model a stokeslet cyclinder flow

%% Constants.

% Space parameters.
Xsi = 2; Xlen = 500; % X dimensionless length, number of points, resp.
Ysi = 2; Ylen = 500; % Y dimensionless length, number of points, resp.
Xsys = linspace(-Xsi/2,Xsi/2,Xlen);
Ysys = linspace(-Ysi/2,Ysi/2,Ylen);

% Store for final velocity field.
UxSys = zeros(Xlen,Ylen);
UySys = zeros(Xlen,Ylen);

% Cylinder parameters
Nstok = 2; % Stokeslets per cylinder.
R = 0.25; % Radius of cylinder.
eps = 2*pi*R/Nstok; % Regularization parameter.
theta0 = pi/4; % Phase shift for cylinder boundary velocity field.
omega = 4/3; % Periodicity modifier for cylinder boundary velocity field.
u0 = 0; % (Note: u0=0 corresponds to no surface velocity) Maximum velocity of cylinder surface.

% Parameterise the cylinder.
theta = linspace(-pi/2,3*pi/2,Nstok+1)'; % Calculate theta to parameterise the cyclinder surface.
theta = theta(1:end-1);

% Get the stokeslet cartesian coordinates.
X = R*cos(theta); % Get the X-coord of the stokeslets.
Y = R*sin(theta); % Get the Y-coord of the stokeslets.

% % Calculate the surface velocities on the cylinder.
% uTmag = zeros(Nstok,1); % Array of magnitudes.
 uTx = zeros(Nstok,1); % X-components of velocity.
 uTy = zeros(Nstok,1); % Y-components of velocity.
% a = find( theta < pi ); thetaA = theta(a); % Get the angles for which the velocity is non-zeros.
% uTmag(a) = u0*sin(omega*(thetaA - theta0)); % Calculate the velocity magnitudes.
% uTx(a) =  -uTmag(a).*sin(thetaA); uTy(a) = uTmag(a).*cos(thetaA); % Calculate the velocity components (tangent to the surface).
 uLen = 2*length(uTx);
 u = zeros(uLen,1);
% 
% % Produce the 'u' column vector to inverse solve for the forces.
% for i = 1:Nstok
%     u(2*i-1) = uTx(i);u(2*i) = uTy(i); % Get the components for the column vector
% end

%% (OPTIONAL) Circle moving in a 2D flow U0 = (0,1).

% Set the velocity so at the surface 
uTx(:) = -1; % No motion in the x direction.
uTy(:) = 0; % Surface velocity cancels the background flow.

u = zeros(uLen,1);

% Produce the 'u' column vector to inverse solve for the forces.
for i = 1:Nstok
    u(2*i-1) = uTx(i);u(2*i) = uTy(i); % Get the components for the column vector
end

%% Calculate the Stokes Matrix.
% At the moment this code is 'additive' and S is designed to change size. A
% more optimal (precallocation) approach would be better, particularly for large Nstok values.

S = []; % Initialise store for final matrix.
Stemp = []; % Initilise store for 2x2 Stokes matrix produced in each comparison.

for i = 1:Nstok % Loop through stokeslets.
    Srowtemp = []; % (Re)initialise the row matrix of (+)_N (2x2).
    for j = 1:Nstok % Loop through stokelets to compare.
        Xstok = [X(i),Y(i)]; % Get Coord of stokeslet 1.
        Ystok = [X(j),Y(j)]; % Get Coord of stokeslet 2.
        Stemp = regStok2D(Xstok,Ystok,eps,1/(4*pi)); % Get the 2x2 associated to X, Y.
        Srowtemp = [Srowtemp, Stemp]; % Build the row of matrix associated to the ith stokeslet.
    end
    S = [S;Srowtemp]; % Add the row to the full matrix.
end

S(S(:,:)<1e-3)=0;

%% Calculate the forces on the boundary.

F = u'/S; % Use the velocity boundary conditions and the stokes matrix to calculate the forces.

% Extract the force components.
Fx = zeros(Nstok,1); % X-components of velocity.
Fy = zeros(Nstok,1); % Y-components of velocity.
for i = 1:Nstok
    Fx(i) = F(2*i-1);
    Fy(i) = F(2*i);
end

% Enforce zero net force
FSX = sum(Fx);
FSY = sum(Fy);
%Fx = Fx - 2*FSX/Nstok;
%Fy = Fy - 2*FSY/Nstok;

%% Calculate the flow field for the whole space.

% Loop over velocity field calculation points.
for xi = 1:Xlen
    for yi = 1:Ylen

        pos = [Xsys(xi),Ysys(yi)]; % Position at which flow is being calculated.
        Utemp = [0,0]'; % Reinitialise the local velocity field storage.

        % Loop over stokeslets.
        for n = 1:Nstok
            posStok = [X(n),Y(n)]; % Position of stokeslet (force) being considered.
            Utemp = Utemp + regStok2D(pos,posStok,eps,1/(4*pi))*[Fx(n),Fy(n)]'; % Stokelet.
        end

        Utemp = Utemp + [1,0]'; % Background flow.

        % Extract local components for velocity field.
        UxSys(xi,yi) = Utemp(1);
        UySys(xi,yi) = Utemp(2);

    end
end

%% plotting tools

%plot(theta,uTmag)
%quiver(X,Y,uTx,uTy)

%% plotting tools 2

contour(Xsys,Ysys,UxSys',100)
%contour(UySys)
axis equal

%% Functions

% Inputs: 
% - x, location of target.
% - y, location of comparison stokelet.
% - eps, regularisation parameter.
% - fpm, 1/4*pi*mu.

function [S] = regStok2D(x,y,eps,fpm)

R = sqrt(norm(x-y) + eps^2) + eps;
rho = (R+eps)/(R*(R-eps));

S = zeros(2,2);

for i = 1:2
    for j = 1:2
       S(i,j) = -isequal(i,j)*(log(R) - eps*rho) + (x(i)-y(i))*(x(j)-y(j))*rho/R;
    end
end

S = S*fpm;

end

% TBD - not functional in current state.
function [] = regStok3D(x,y,eps)

R = sqrt(norm(x-y) + eps^2);

S = zeros(3,3);

for i = 1:3
    for j = 1:3
        S(i,j) = isequal(i,j)*(1/R + eps^2/R^3) + (x(i)-y(i))*(x(j)-y(j))/R^3;
    end
end

end