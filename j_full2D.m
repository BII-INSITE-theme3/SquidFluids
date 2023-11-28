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




%% Notes, based on old code, remove when working

% hold on;
% c = jet(max(stks(:,3))); % (Optional)
% scatter(stks(:,1),stks(:,2),2,c(stks(:,3),:)) % (Optional)
% axis equal % (Optional)

%%

% Put the fluid velocity into "vertical array" form.
for i = 1:Nstoks
    U0V(2*i-1) = U0(i,2);
    U0V(2*i)   = U0(i,1);
end

% Find the forces.
for l = 1:2*Nstoks % Loop through the stokeslet-components
    for k = 1:2*Nstoks % Loop through the stokeslet-components
        %
        n1 = ceil(k/2); % Get the stokelets number of the influence stokeslet.
        n2 = ceil(l/2); % Get the stokelets number of the influenced stokeslet.
        %
        stks1 = stks(n1,:); % Get the position of the influence stokeslet.
        stks2 = stks(n2,:); % Get the position of the influenced stokeslet.
        %
        r = sqrt(norm(stks1-stks2)^2 + eps^2) + eps; % Get the reg `distance' between them.
        rho = (r+eps)/(r*(r-eps)); % Get the "rho", for convenience.
        %
        S(l,k) = -(log(r)-eps*rho)*(mod(k,2)==mod(l,2)) ... % Get the log term contribution.
            + (stks1(1+mod(k,2))-stks2(1+mod(k,2)))*(stks1(1+mod(l,2))-stks2(1+mod(l,2)))*rho/r; % Get the <....> contribution.
    end
end

FV = U0V/S; % Get the forces from the stokeslet.

% Put the forces velocity into vertical array form.
for i = 1:Nstoks
    F(i,2) = FV(2*i-1);
    F(i,1) = FV(2*i);
end

% Run to here for the forces.

%% Solve over the space.

% Loop over the whole space.
parfor i = 1:Npts
    for j = 1:Npts

        Stemp = zeros(2,2);
        tempStks = stks;
        tempF = F;
        tempX = x;
        tempY = y;

        p = [tempX(j),tempY(i)]; % Get the position of consideration, i and j done make sense to me here??

        for n = 1:Nstoks

            pN = tempStks(n,:); % Get the position of stokeslet N.
            Ftemp = tempF(n,:); % Get the forces of stokeslet N.
            r = sqrt(norm(p - pN).^2 + eps^2) + eps; % Distance, considered to stokeslet N.
            %r = sqrt(norm(p - pN).^2 + eps^2) + eps; % Distance, considered to stokeslet N.
            rho = (r+eps)/(r*(r-eps)); % Rho, considered to stokeslet N.

            for k = 1:2
                for l = 1:2
                    Stemp(k,l) = -(log(r)-eps*rho)*(k==l) + (p(k)-stks(n,k))*(p(l)-stks(n,l))*rho/r;
                    %Stemp(k,l) = -(log(r)-eps*rho)*(k==l) + (p(k)-pN(k))*(p(l)-pN(l))*rho/r;
                end
            end

            U = Stemp*Ftemp';
            Ux(i,j) = Ux(i,j) + U(1);
            Uy(i,j) = Uy(i,j) + U(2);

        end

    end
end

%%

n = 5;

x = linspace(-Xma/2,Xma/2,Npts); % X-values of solution space.
y = linspace(-Yma/2,Yma/2,Npts); % Y-values of solution space.

Uxtemp = Ux + Uflow(1);
Uytemp = Uy + Uflow(2);

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
