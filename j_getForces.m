% Function to solve the inverse problem to get the stokeslet forces.

function [F] = j_getForces(stks)

    F=1; % To delete

    % Create the array of boundary velocities (from BVP)
    [nStok,~] = size(stks);
    BdryVelo = zeros(nStok,2);

    %% No-slip boundaries -- stks(:,3) == 1

    ind = find(stks(:,3)==1);
    BdryVelo(ind,:) = 0; % Set zero-flow on the channel boundaries

    %% No-slip boundaries -- stks(:,3) == 2

    ind = find(stks(:,3)==2);
    nTemp = length(ind);
    BdryVelo(ind,:) = j_poisuelleFlow(nTemp);

    %% No-slip boundaries -- stks(:,3) == 3
    % Can add other boundary conditions for the bottom of a channel here

    %% No-slip boundaries -- stks(:,3) == 4

    ind = find(stks(:,3)==4);nTemp = length(ind);
    BdryVelo(ind,:) = j_surfaceFlow(nTemp,floor(0.75*nTemp));

    ind = find(stks(:,3)==5);nTemp = length(ind);
    BdryVelo(ind,:) = j_surfaceFlow(nTemp,floor(0.5*nTemp));

    ind = find(stks(:,3)==6);nTemp = length(ind);
    BdryVelo(ind,:) = j_surfaceFlow(nTemp,floor(0.5*nTemp));

    ind = find(stks(:,3)==7);nTemp = length(ind);
    BdryVelo(ind,:) = j_surfaceFlow(nTemp,floor(0.75*nTemp));

    % quiver(stks(:,1),stks(:,2),BdryVelo(:,1),BdryVelo(:,2))
    % axis equal

end

