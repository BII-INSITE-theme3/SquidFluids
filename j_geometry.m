% Function to initialise the geometry of the system.
% Outputs an Nx2 array with N (x,y) positisons for each of the stokeslets.

function stks = j_geometry(rho1,Lt1,Lm1,Lb1,theta1,Ptx1,Pty1,dsep1,psi1,PRAx1,PRAy1)

    stks_channel = j_geometry_poisuelle(rho1,Lt1,Lm1,Lb1,theta1,Ptx1,Pty1); % Set the channel geometry.

    %stks_cap = j_geometry_capsule(); % Set the channel geometry with a capsule.

    stks_appendages1 = j_geometry_cylinderPair(rho1,dsep1,psi1,PRAx1,PRAy1,0); % Set the left appendage pair geometry.

    stks_appendages2 = j_geometry_cylinderPair(rho1,dsep1,pi-psi1,-PRAx1,PRAy1,1); % Set the left appendage pair geometry.

    stks = [stks_channel;stks_appendages1; stks_appendages2]; % Combine all structures.
    %stks = [stks_channel;stls_cap; stks_appendages1; stks_appendages2]; % Combine all structures.
    
end