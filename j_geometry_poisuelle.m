% Function to initialise the geometry of the channel with a contraction.
% Outputs an Nx3 array with N (x,y) positisons for each of the stokeslets, with the third argument a one-hot encoded boundary condition.
% boundary codes: (:,3) = 1, no-slip/pen. (:,3) = 2, Poisuelle. (:,4) = 3, free.

function [channel_stks] = j_geometry_poisuelle(rho2,Lt2,Lm2,Lb2,theta2,Ptx2,Pty2)
    
    removeNum = 1;

    % Construct the right boundary.
    WRTl = linspace(0,Lt2,floor(Lt2*rho2)+removeNum); % Parameterise the right wall top section.
    WRTl = WRTl(1:end-removeNum);
    WRMl = linspace(0,Lm2,floor(Lm2*rho2)); % Parameterise the right wall middle section.
    WRBl = linspace(0,Lb2,floor(Lb2*rho2)+removeNum); % Parameterise the right wall bottom section.
    WRBl = WRBl(1+removeNum:end);

    WRT = [Ptx2 - WRTl*sin(0); Pty2 - WRTl*cos(0);ones(floor(Lt2*rho2),1)']'; % Get the top segment coordinates.
    WRM = [Ptx2 - WRMl*sin(theta2); Pty2 - WRTl(end)*cos(0) - WRMl*cos(theta2);ones(floor(Lm2*rho2),1)']'; % Get the middle segment coordinates.
    WRB = [Ptx2 - WRMl(end)*sin(theta2) - WRBl*sin(0);Pty2 - WRTl(end)*cos(0) - WRMl(end)*cos(theta2) - WRBl*cos(0);ones(floor(Lb2*rho2),1)']'; % Get the bottom segment coordinates.
    
    wallR = [WRT;WRM;WRB]; % Combine all the segments to give the right wall.

    % Construct the left boundary.
    wallL = wallR; % Copy the left boundary coordinates.
    wallL(:,1) = -wallL(:,1); % Reflect x-coordinate.

    % Construct the top boundary.
    wallT = [linspace(-Ptx2,Ptx2,floor(2*Ptx2*rho2)); Pty2*ones(floor(2*Ptx2*rho2),1)'; 2*ones(floor(2*Ptx2*rho2),1)']';
    
    % Construct the bottom boundary.
    % Lbottom = Ptx2 - WRMl(end)*sin(theta2); % x-distance from horizontal to bottom segement.
    % Ybottom = Pty2 - WRTl(end)*cos(0) - WRMl(end)*cos(theta2) - WRBl(end)*cos(0);  % y-coordinate to bottom segement.
    % wallB = [linspace(-Lbottom,Lbottom,floor(2*Lbottom*rho2)); Ybottom*ones(floor(2*Lbottom*rho2),1)'; 3*ones(floor(2*Lbottom*rho2),1)']';

    % Combine to give the full position array.
    channel_stks = [wallL;wallR;wallT];
    %channel_stks = [wallL;wallR;wallT;wallB];

end

