% Title: Second attempt to do the inverse problem.
% Author: Stephen Williams.
% Notes: ...
%--------------------------------------------%

% Constants.
Xma = 1;
Yma = Xma;
Npts = 50; % Number of points in the "calculated" space.
eps = 0.25*Xma/(Npts); % Epsilon of the regularization.
R0 = 0.25; % Radius of the stokelet surface.
Nstoks = 160; % Number of stokeslets on the radius.

% Preallocation.
x = linspace(-Xma/2,Xma/2,Npts); % X-values of solution space.
y = linspace(-Yma/2,Yma/2,Npts); % Y-values of solution space.
%
Ux = zeros(Npts); % fluid velocity x-component.
Uy = zeros(Npts); % fluid velocity y-component.
%
theta = linspace(0,2*pi,Nstoks+1); % Get the angles on the surface of the stokeslets.
theta = theta(1:end-1); % Remove the repeat value.
%
stks = R0*[cos(theta'),sin(theta');]; % Get the stokeslet cartesian coordinates.
%
S = zeros(2*Nstoks); % Store for the Stokeslet.
Stemp = zeros(2); % Store for the Stokeslet.
F = zeros(Nstoks,2); % Store for the stokeslet's forces.
Ftemp = zeros(2,1); % Store for the stokeslet's forces.
%
U0 = zeros(Nstoks,2); % Store the for fluid velocity at the stokeslets.
U0V = zeros(1,2*Nstoks); % Vertical store the for fluid velocity at the stokeslets.

% Boundary flow field.
U0(:,1) = 1; % Set the values (this is where the code is modified).

% Put the fluid velocity into "vertical array" form.
for i = 1:Nstoks
    U0V(2*i-1) = U0(i,1);
    U0V(2*i)   = U0(i,2);
end

% Find the forces.
for k = 1:2*Nstoks % Loop through the stokeslet-components
    for l = 1:2*Nstoks % Loop through the stokeslet-components
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
        S(k,l) = -(log(r)-eps*rho)*(mod(k,2)==mod(l,2)) ... % Get the log term contribution.
        + (stks1(1+mod(k,2))-stks2(1+mod(k,2)))*(stks1(1+mod(l,2))-stks2(1+mod(l,2)))*rho/r; % Get the <....> contribution.
    end
end

FV = U0V/S; % Get the forces from the stokeslet.

% Put the forces velocity into vertical array form.
for i = 1:Nstoks
    F(i,1) = FV(2*i-1);
    F(i,2) = FV(2*i);
end

%% Solve over the space.

% Loop over the whole space.
for i = 1:Npts
    for j = 1:Npts

        p = [x(i),y(j)]; % Get the position of consideration

        for n = 1:Nstoks

            pN = stks(n,:); % Get the position of stokeslet N.
            Ftemp = F(n,:); % Get the forces of stokeslet N.
            r = sqrt(norm(p' - pN').^2 + eps^2) + eps; % Distance, considered to stokeslet N.
            rho = (r+eps)/(r*(r-eps)); % Rho, considered to stokeslet N.

            for k = 1:2
                for l = 1:2
                    Stemp(k,l) = -(log(r)-eps*rho)*(k==l) + (p(k)-stks(n,k))*(p(l)-stks(n,l))*rho/r;
                end
            end

            U = Stemp*Ftemp';

            Ux(i,j) = Ux(i,j) + U(1);
            Uy(i,j) = Uy(i,j) + U(2);

        end
    end
end

%%

n = 100;

UxTemp = Ux;
UyTemp = Uy;

% Optional thresholding to improve the contour plots.
thresh1 = 10;
UxTemp( abs(Ux) > thresh1) = thresh1;
UyTemp( abs(Uy) > thresh1) = thresh1;

figure
contour(Xsys,Ysys,UySys',n,'k')
hold on
contour(x,y,UyTemp',n/10,'r')

%quiver(x(1:n:end),y(1:n:end),Ux(1:n:end,1:n:end),Uy(1:n:end,1:n:end),5)
