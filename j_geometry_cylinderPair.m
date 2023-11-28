% Function to initialise the geometry of the appendages.
% boundary codes: (:,3) = 4, left of right pair. boundary codes: (:,3) = 5, right of right pair.

function [stks_cyl] = j_geometry_cylinderPair(rho2,dsep2,psi2,PRAx2,PRAy2,mod1)
    
    theta = linspace(0,2*pi,floor(rho2*2*pi)+1); % Get the angles to parameterise the surfaces.
    theta = theta(1:end-1); % Remove repeat.
    
    circle = [sin(theta);cos(theta);zeros(floor(rho2*2*pi),1)']';

    stks_cyl1 = circle + [PRAx2+(1+dsep2/2)*cos(psi2),PRAy2+(1+dsep2/2)*sin(psi2),4+2*mod1];
    stks_cyl2 = circle + [PRAx2-(1+dsep2/2)*cos(psi2),PRAy2-(1+dsep2/2)*sin(psi2),5+2*mod1];

    stks_cyl  = [stks_cyl1;stks_cyl2];

end

