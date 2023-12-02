% Title: Second attempt to do the inverse problem.
% Author: Stephen Williams.
% Notes: 1. Code is non-dimensionalised. 
% 1.1 The characteristic length is equal to the appendage radius. 
% 1.2 The characteristic time is the time it would take for the flow at the maximum on the entry point to move 1 appendage radius.
% 1.3 The viscocity is set so that the prefactor in the stokeslet is 1.
% 2. The boundaries are now arc-length parameterisation (instead of uniform
% along y). This improves the stokeslet distribution on the non-vertical
% portions of the boundaries.
%--------------------------------------------%

close all
clear all

%% Set parameters
j_parameters % Set the parameters, TBD

%% Initialise parallelisation
% parpool % Intialise the parallel workers % (Optional)

%% Set channel geometry
% Set the system geometry
stks = j_geometry(rho,Lt,Lm,Lb,theta,Ptx,Pty,dsep,psi,PRAx,PRAy);

%% Solve for the forces
F = j_getForces(stks,eps,U0);

%% Check the boundary flow is as expected

n=1;
[Uflowx,Uflowy] = j_calculateFlowPath(stks,F,stks(:,1),stks(:,2));
quiver(stks(:,1),stks(:,2),Uflowx(1:n:end),Uflowy(1:n:end)) % Plot the vector field

%% Get the flow over the space
 [Uflowx,Uflowy] = j_calculateFlowGrid(stks,F,x,y);

%% Plot output flow field

Uflowy = Uflowy - mean(mean(Uflowy));
Uflowx = Uflowx - mean(mean(Uflowx));

n=1; % Plot coarseness
Umag = sqrt(Uflowx.^2 + Uflowy.^2); % Get flow field magnitude
imagesc(x,y,Umag') % Plot local flow strenght as a background
hold on
scatter(stks(:,1),stks(:,2),2,'r') % Plot the Stokeslets
quiver(x(1:n:end),y(1:n:end),Uflowx(1:n:end,1:n:end),Uflowy(1:n:end,1:n:end)) % Plot the vector field
axis equal
