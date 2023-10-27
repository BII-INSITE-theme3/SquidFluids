
% Title: 2D cylinder pair code (refactor).
% Author: Stephen Williams.
% Notes: 
%--------------------------------------------%

%% Constants and preallocation

% System geometry
Npts = 100; % Number of points in the space to calculate.
Xma = 1; % Largest x to calculate (system symmetric).
Yma = Xma; % Largest y to calculate (system symmetric).
x = linspace(-Xma,Xma,Npts); % x coordinates of the fluid calculations.
y = linspace(-Yma,Yma,Npts); % y coordinates of the fluid calculations.

% Stokeslet parameters
eps = 1e-5; % Stokeslet regularisation parameter.
Nstructures = 2; % Number of structures to construct.
Nstokeslets = 80*ones(Nstructures,1); % Number of Stokelsets per structure.
NStokesTotal = sum(Nstokeslets);
structCenters = [-0.3,0; 0.3,0]; % Positions of the structure centers.
structRadii = [0.2,0.2]; % Radii of the structures.
PosStokeslets = setPositions(Nstructures,Nstokeslets,structCenters,structRadii);
S = zeros(2*NStokesTotal); % Store for the overall strucutre Stokeslet.
Stemp = zeros(2); % Store for temporary Stokeslets.
F = zeros(NStokesTotal,2);
Ftemp = zeros(2,1);

% Fluid parameters
Ux = zeros(Npts);Uy = zeros(Npts); % Stores for the fluid velocities.
bckgrdang = pi/2; % Background flow angle to horizontal, anticlockwise.
bckgrdstr = 1; % Background flow strength.
Ubackgrd = bckgrdstr*[cos(bckgrdang),sin(bckgrdang)]; % Background flow.

%% Set fluid flow on the structure boundaries.

Ubdry = setBoundaryFlow(Ubackgrd,Nstructures,Nstokeslets);

%% Calculate the forces



%% Use the forces to calculate the flow



%% Plotting tools



%% Functions

% Get the positions of all of the Stokeslets in the system
% To do this, we assume  they are uniformly distributed on the surfaces of
% the structures, which are themselves assumed to be circular.
function [pos] = setPositions(Nstructures,Nstokeslets,structCenters,structRadii)

    pos = []; % Array of positions to be output.

    % Loop over the structures.
    for i = 1:Nstructures

        theta = linspace(0,2*pi,Nstokeslets(i)+1); % Get the uniformly distributed angles.
        theta = theta(1:end-1); % Remove the repeat.
        posi = [structCenters(i,1) + structRadii(i)*cos(theta); ... % Get the x positions.
                structCenters(i,2) + structRadii(i)*sin(theta)]; % Get the y positions.
        pos = [pos,posi]; % Append to position array

    end
end

function [Ubdry] = setBoundaryFlow(Ubackgrd,Nstructures,Nstokeslets)

    for i = 1:Nstructures

        

    end

end
