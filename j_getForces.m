% Function to solve the inverse problem to get the stokeslet forces.

function [F] = j_getForces(stks1,eps1,U01)

    % Create the array of boundary velocities (from BVP)
    [nStok,~] = size(stks1);
    BdryVelo = zeros(nStok,2);

    %% No-slip boundaries -- stks(:,3) == 1, this code currently does nothing, so is commented out.
    %ind = find(stks(:,3)==1);
    %BdryVelo(ind,:) = 0; % Set zero-flow on the channel boundaries

    %% Poisuelle boundary flow -- stks(:,3) == 2
    % Flow going as 1-(r/a)^2, r = distance from channel center, a = half channel-width.

    ind = find(stks1(:,3)==2);
    nTemp = length(ind); 
    BdryVelo(ind,:) = j_poisuelleFlow(nTemp,U01);

    %% No-slip boundaries -- stks(:,3) == 3
    % Can add other boundary conditions for the bottom of a channel here.

    %% No-slip boundaries -- stks(:,3) == 4,5,6,7
    % Boundary conditions following those in Nawroth 2017.

    ind = find(stks1(:,3)==4); nTemp = length(ind);
    BdryVelo(ind,:) = j_surfaceFlow(nTemp,floor(0.75*nTemp));

    ind = find(stks1(:,3)==5); nTemp = length(ind);
    BdryVelo(ind,:) = j_surfaceFlow(nTemp,floor(0.5*nTemp));

    ind = find(stks1(:,3)==6); nTemp = length(ind);
    BdryVelo(ind,:) = j_surfaceFlow(nTemp,floor(0.5*nTemp));

    ind = find(stks1(:,3)==7); nTemp = length(ind);
    BdryVelo(ind,:) = j_surfaceFlow(nTemp,floor(0.75*nTemp));

    %% Create the linear system to solve, part 1.
    % Put the velocity boundary conditions in vertical form.

    BdryVertical = zeros(nStok*2,1);
    BdryVertical(1:2:end) = BdryVelo(:,1); BdryVertical(2:2:end) = BdryVelo(:,2);
    % quiver(stks1(:,1),stks1(:,2),BdryVelo(:,1),BdryVelo(:,2))

    %% Create the linear system to solve, part 2.
    % Create a matrix corresponding to the stokeslet for the full system.

    S = zeros(2*nStok); % Preallocate the full linear stokeslet.

    % Find the forces.
    for l = 1:2*nStok % Loop through the stokeslet-components
        for k = 1:2*nStok % Loop through the stokeslet-components
            %
            n1 = ceil(k/2); % Get the stokelets number of the influence stokeslet.
            n2 = ceil(l/2); % Get the stokelets number of the influenced stokeslet.
            %
            tempStks1 = stks1(n1,:); % Get the position of the influence stokeslet.
            tempStks2 = stks1(n2,:); % Get the position of the influenced stokeslet.
            %
            r = sqrt(norm(tempStks1-tempStks2)^2 + eps1^2) + eps1; % Get the reg `distance' between them.
            rho = (r+eps1)/(r*(r-eps1)); % Get the "rho", for convenience.
            %
            S(l,k) = -(log(r)-eps1*rho)*(mod(k,2)==mod(l,2)) ... % Get the log term contribution.
                + (tempStks1(1+mod(k,2))-tempStks2(1+mod(k,2)))*(tempStks1(1+mod(l,2))-tempStks2(1+mod(l,2)))*rho/r; % Get the <....> contribution.
        end
    end

    %% Create the linear system to solve, part 2.
    % Solve for the forces required to satisfy the system.

    % ForceVertical = zeros(2*nStok,1);
    ForceVertical = (BdryVertical\S);
    F = zeros(nStok,2);
    F(:,1) = ForceVertical(1:2:end);
    F(:,2) = ForceVertical(2:2:end);

    % F(abs(F(:,1))>4,1)=4; F(abs(F(:,2))>4,2)=4;
    quiver(stks1(:,1),stks1(:,2),F(:,1),F(:,2))

end

