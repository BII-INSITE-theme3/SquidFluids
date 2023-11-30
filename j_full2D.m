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
% Set the system geometry.
stks = j_geometry(rho,Lt,Lm,Lb,theta,Ptx,Pty,dsep,psi,PRAx,PRAy);

%% Solve for the forces
F = j_getForces(stks);

%% Get the flow over the space
[Uflowx,Uflowy] = j_calculateFlow(stks,F,x,y);

%% Notes, based on old code, remove when working

% hold on;
% c = jet(max(stks(:,3))); % (Optional)
% scatter(stks(:,1),stks(:,2),2,c(stks(:,3),:)) % (Optional)
% axis equal % (Optional)

%%

Uxtemp = Uflowx;
Uytemp = Uflowy;

% Optional thresholding to improve the contour plots.
% thresh1 = 20;
% UxTemp( abs(Ux) > thresh1) = thresh1;
% UyTemp( abs(Uy) > thresh1) = thresh1;

Umag = sqrt(Uxtemp.^2 + Uytemp.^2);

figure
imagesc(x,y,Umag)
%set(gca,'YDir','normal')
hold on
%contour(x,y,UxTemp',n,'r')
scatter(stks(:,1),stks(:,2),2,'r')
quiver(x(1:n:end),y(1:n:end),Uxtemp(1:n:end,1:n:end),Uytemp(1:n:end,1:n:end),2)
axis equal
