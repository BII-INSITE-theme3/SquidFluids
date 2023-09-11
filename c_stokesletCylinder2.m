% Title: Second attempt to do the inverse problem.
% Author: Stephen Williams.
% Notes: ...
%--------------------------------------------%

% Constants.
Npts = 100; % Number of points in the "calculated" space.
eps = 0.0001; % Epsilon of the regularization.
R0 = 0.2; % Radius of the stokelet surface.
Nstoks = 5; % Number of stokeslets on the radius.

% Preallocation.
x = linspace(-1,1,Npts); % X-values of solution space.
y = linspace(-1,1,Npts); % Y-values of solution space.
Ux = zeros(Npts); % fluid velocity x-component.
Uy = zeros(Npts); % fluid velocity y-component.
theta = linspace(0,2*pi,Nstoks+1); % Get the angles on the surface of the stokeslets.
theta = theta(1:end-1); % Remove the repeat value.
stks = R0*[cos(theta'),sin(theta');]; % Get the cartesian coordinates.
S = zeros(2*Nstoks); % Store for the Stokeslet.
U0 = zeros(Nstoks,2); % Store the for fluid velocity at the stokeslets.
U0V = zeros(1,2*Nstoks); % Vertical store the for fluid velocity at the stokeslets.
U0(:,1) = 1; % Set the values (this is where the code is modified).

% Put the fluid velocity into "vertical array" form.
for i = 1:Nstoks
    U0V(2*i-1) = U0(i,1);
    U0V(2*i)   = U0(i,2);
end

% Find the forces.
for k = 1:2*Nstoks % Loop through the stokeslet-components
    for l = 1:2*Nstoks % Loop through the stokeslet-components

        n1 = ceil(k/2); % Get the stokelets number of the influence stokeslet.
        n2 = ceil(l/2); % Get the stokelets number of the influenced stokeslet.

        stks1 = stks(n1,:); % Get the position of the influence stokeslet.
        stks2 = stks(n2,:); % Get the position of the influenced stokeslet.

        r = sqrt(norm(stks1-stks2) + eps^2) + eps; % Get the distance between them.
        rho = (r+eps)/(r*(r-eps)); % Get the "rho", for convenience.

        S(k,l) = -(log(r)-eps*rho)*(mod(k,2)==mod(l,2)) ... % Get the log term contribution.
        + (stks1(1+mod(k,2))-stks2(1+mod(k,2)))*(stks1(1+mod(l,2))-stks2(1+mod(l,2)))*rho/r; % Get the <....> contribution.
    end
end

FV = U0V/S; % Get the forces from the stokeslet.

% Put the forces velocity into vertical array form.
F = zeros(Nstoks,2);
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

            

        end
        
        % for k = 1:2*Nstok
        %     for l = 1:2*Nstok
        % 
        %         stokN = 123;
        % 
        %         r = sqrt(norm(p,stks) + eps^2) + eps;
        %         rho = R+eps/(R*(R-eps));
        % 
        %         S(k,l) = (log(r)-eps*rho)*(k==l) + (p(k)-stks(k))*(p(j)-stks(j))*rho/r;
        %     end
        % end
        
        U = S*F; % Get the fluid velocity at the considered point.
        Ux(i,j) = U(1); % Get the fluid velocity x-component at the considered point.
        Uy(i,j) = U(2); % Get the fluid velocity y-component at the considered point.

    end
end
